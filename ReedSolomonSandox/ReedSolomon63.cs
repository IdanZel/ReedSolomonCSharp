namespace ReedSolomonSandox;

public class ReedSolomon63 : ReedSolomonTabs
{
    private const int NROOTS = 27;
    private const uint LOAD = SIZE - NROOTS;
    private const int FCR = 1;

    public unsafe int Decode(byte[] data, uint len, byte[]? parity, uint[] erasPos, uint noEras, byte[] corr)
    {
        if (len < (parity is not null ? 1 : NROOTS + 1))
        {
            throw new Exception("reed-solomon: must provide all parity and at least one non-parity symbol");
        }

        if (parity is null)
        {
            len -= NROOTS;
            parity = data.Skip((int)len).ToArray();
        }

        byte[] tmpArr = new byte[SIZE];

        fixed (byte* tmp = tmpArr)
        {
            byte mask;
            unchecked
            {
                mask = (byte)(~0L << SYM);
            }

            int dataStartIndex = (int)(LOAD - len);

            for (int i = 0; i < len; i++)
            {
                tmp[dataStartIndex + i] = (byte)(data[i] & ~mask);
            }

            byte* dataBytes = tmp + dataStartIndex;

            for (int i = 0; i < NROOTS; i++)
            {
                if ((parity[i] & mask) != 0)
                {
                    throw new Exception("reed-solomon: parity data contains information beyond R-S symbol size");
                }

                tmp[LOAD + i] = parity[i];
            }

            byte* parityBytes = tmp + LOAD;

            int corrects = DecodeSymbols(dataBytes, len, parityBytes, erasPos, noEras, corr);

            if (corrects > 0)
            {
                for (int i = 0; i < len; i++)
                {
                    data[i] &= mask;
                    data[i] |= tmp[LOAD - len + i];
                }

                for (int i = 0; i < NROOTS; i++)
                {
                    parity[i] &= tmp[LOAD + i];
                }
            }

            return corrects;
        }
    }

    private unsafe int DecodeSymbols(
        byte* data,
        uint len,
        byte* parity, // Requires: at least NROOTS
        uint[] erasPos, // Capacity: at least NROOTS
        uint noEras, // Maximum:  at most  NROOTS
        byte[] corr) // Capacity: at least NROOTS
    {
        byte[] lambda = new byte[NROOTS + 1];
        byte[] syn = new byte[NROOTS];
        byte[] b = new byte[NROOTS + 1];
        byte[] t = new byte[NROOTS + 1];
        byte[] omega = new byte[NROOTS + 1];
        uint[] root = new uint[NROOTS];
        uint[] loc = new uint[NROOTS];
        int count;

        // Check length parameter for validity.  We know NROOTS < NN.  It is possible to have as
        // little as 1 non-parity (payload) symbol that isn't a pad, and as many as LOAD (SIZE -
        // NROOTS.
        if (len == 0 || len > LOAD)
        {
            throw new Exception("reed-solomon: data length incompatible with block size and error correction symbols");
        }

        uint pad = LOAD - len;
        if (noEras > 0)
        {
            if (noEras > NROOTS)
            {
                throw new Exception("reed-solomon: number of erasures exceeds capacity (number of roots)");
            }

            for (int i = 0; i < noEras; ++i)
            {
                if (erasPos[i] >= len + NROOTS)
                {
                    throw new Exception("reed-solomon: erasure positions outside data+parity");
                }
            }
        }

        // form the syndromes; i.e., evaluate data(x) at roots of g(x)
        for (int i = 0; i < NROOTS; i++)
            syn[i] = data[0];

        for (int j = 1; j < len; j++)
        {
            for (int i = 0; i < NROOTS; i++)
            {
                if (syn[i] == 0)
                {
                    syn[i] = data[j];
                }
                else
                {
                    syn[i] = (byte)(data[j] ^ AlphaTo[Modnn((uint)(IndexOf[syn[i]] + (FCR + i) * PRM))]);
                }
            }
        }

        for (int j = 0; j < NROOTS; j++)
        {
            for (int i = 0; i < NROOTS; i++)
            {
                if (syn[i] == 0)
                {
                    syn[i] = parity[j];
                }
                else
                {
                    syn[i] = (byte)(parity[j] ^ AlphaTo[Modnn((uint)(IndexOf[syn[i]] + (FCR + i) * PRM))]);
                }
            }
        }

        // Convert syndromes to index form, checking for nonzero condition
        byte synError = 0;
        for (int i = 0; i < NROOTS; i++)
        {
            synError |= syn[i];
            syn[i] = IndexOf[syn[i]];
        }

        uint degLambda = 0;
        uint r = noEras;
        uint el = noEras;
        if (synError == 0)
        {
            // if syndrome is zero, data[] is a codeword and there are no errors to correct.

            count = 0;
            goto finish;
        }

        lambda[0] = 1;
        if (noEras > 0)
        {
            // Init lambda to be the erasure locator polynomial.  Convert erasure positions
            // from index into data, to index into Reed-Solomon block.
            lambda[1] = AlphaTo[Modnn(PRM * (NN - 1 - (erasPos[0] + pad)))];
            for (int i = 1; i < noEras; i++)
            {
                byte u = Modnn(PRM * (NN - 1 - (erasPos[i] + pad)));
                for (int j = i + 1; j > 0; j--)
                {
                    byte tmp = IndexOf[lambda[j - 1]];
                    if (tmp != A0)
                    {
                        lambda[j] ^= AlphaTo[Modnn((uint)(u + tmp))];
                    }
                }
            }
        }

        for (int i = 0; i < NROOTS + 1; i++)
            b[i] = IndexOf[lambda[i]];

        //
        // Begin Berlekamp-Massey algorithm to determine error+erasure locator polynomial
        //
        while (++r <= NROOTS)
        {
            // r is the step number
            // Compute discrepancy at the r-th step in poly-form
            byte discrR = 0;
            for (int i = 0; i < r; i++)
            {
                if ((lambda[i] != 0) && (syn[r - i - 1] != A0))
                {
                    discrR ^= AlphaTo[Modnn((uint)(IndexOf[lambda[i]] + syn[r - i - 1]))];
                }
            }

            discrR = IndexOf[discrR]; // Index form
            if (discrR == A0)
            {
                // 2 lines below: B(x) <-- x*B(x)
                // Rotate the last element of b[NROOTS+1] to b[0]
                b = b.Skip(NROOTS).Concat(b.Take(NROOTS)).ToArray();
                b[0] = (byte)A0;
            }
            else
            {
                // 7 lines below: T(x) <-- lambda(x)-discr_r*x*b(x)
                t[0] = lambda[0];
                for (int i = 0; i < NROOTS; i++)
                {
                    if (b[i] != A0)
                    {
                        t[i + 1] = (byte)(lambda[i + 1] ^ AlphaTo[Modnn((uint)(discrR + b[i]))]);
                    }
                    else
                        t[i + 1] = lambda[i + 1];
                }

                if (2 * el <= r + noEras - 1)
                {
                    el = r + noEras - el;
                    // 2 lines below: B(x) <-- inv(discr_r) * lambda(x)
                    for (int i = 0; i <= NROOTS; i++)
                    {
                        b[i] = (byte)(lambda[i] == 0 ? A0 : Modnn((uint)(IndexOf[lambda[i]] - discrR + NN)));
                    }
                }
                else
                {
                    // 2 lines below: B(x) <-- x*B(x)
                    b = b.Skip(NROOTS).Concat(b.Take(NROOTS)).ToArray();
                    b[0] = (byte)A0;
                }

                lambda = t.ToArray();
            }
        }

        // Convert lambda to index form and compute deg(lambda(x))
        for (uint i = 0; i < NROOTS + 1; i++)
        {
            lambda[i] = IndexOf[lambda[i]];
            if (lambda[i] != NN)
                degLambda = i;
        }

        // Find roots of error+erasure locator polynomial by Chien search
        byte[] reg = lambda.ToArray();
        count = 0; // Number of roots of lambda(x)
        for (int i = 1, k = (int)(Iprim - 1); i <= NN; i++, k = Modnn((uint)(k + Iprim)))
        {
            byte q = 1; // lambda[0] is always 0
            for (uint j = degLambda; j > 0; j--)
            {
                if (reg[j] != A0)
                {
                    reg[j] = Modnn(reg[j] + j);
                    q ^= AlphaTo[reg[j]];
                }
            }

            if (q != 0)
                continue; // Not a root
            // store root (index-form) and error location number
            root[count] = (uint)i;
            loc[count] = (uint)k;
            // If we've already found max possible roots, abort the search to save time
            if (++count == (int)degLambda)
                break;
        }

        if ((int)degLambda != count)
        {
            // deg(lambda) unequal to number of roots => uncorrectable error detected
            count = -1;
            goto finish;
        }

        //
        // Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo x**NROOTS). in
        // index form. Also find deg(omega).
        //
        if (degLambda == 0)
        {
            count = -1;
            goto finish;
        }

        var degOmega = degLambda - 1;
        for (int i = 0; i <= degOmega; i++)
        {
            byte tmp = 0;
            for (int j = i + 1; j-- > 0;)
            {
                // int j descending from i to 0, inclusive
                if ((syn[i - j] != A0) && (lambda[j] != A0))
                    tmp ^= AlphaTo[Modnn((uint)(syn[i - j] + lambda[j]))];
            }

            omega[i] = IndexOf[tmp];
        }

        //
        // Compute error values in poly-form. num1 = omega(inv(X(l))), num2 = inv(X(l))**(fcr-1)
        // and den = lambda_pr(inv(X(l))) all in poly-form
        //
        for (int j = count; j-- > 0;)
        {
            byte num1 = 0;
            for (uint i = degOmega + 1; i-- > 0;)
            {
                if (omega[i] != A0)
                    num1 ^= AlphaTo[Modnn(omega[i] + i * root[j])];
            }

            byte num2 = AlphaTo[Modnn(root[j] * (FCR - 1) + NN)];
            byte den = 0;

            // lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i]
            for (int i = (int)(Math.Min(degLambda, NROOTS - 1) & ~1); i >= 0; i -= 2)
            {
                if (lambda[i + 1] != A0)
                {
                    den ^= AlphaTo[Modnn((uint)(lambda[i + 1] + i * root[j]))];
                }
            }

            if (den == 0)
            {
                count = -1;
                goto finish;
            }

            // Apply error to data.  Padding ('pad' unused symbols) begin at index 0.
            if (num1 != 0)
            {
                if (loc[j] < pad)
                {
                    // If the computed error position is in the 'pad' (the unused portion of the
                    // R-S data capacity), then our solution has failed -- we've computed a
                    // correction location outside of the data and parity we've been provided!

                    count = -1;
                    goto finish;
                }

                byte cor = AlphaTo[Modnn((uint)(IndexOf[num1]
                    + IndexOf[num2]
                    + NN - IndexOf[den]))];
                // Store the error correction pattern, if a correction buffer is available.
                // This must be the error correction in the basis of the data/parity buffers,
                // which might be dual-basis encoded.  If so -- convert to conventional, and
                // compute the difference between the erroneous conventional symbol and the
                // corrected conventional.
                if (corr is not null)
                    corr[j] = cor;
                // If a data/parity buffer is given and the error is inside the message or
                // parity data, correct it, converting from and back into dual-basis, if
                // necessary.  The correction will not be in the pad (eliminated, above).
                if (loc[j] < (NN - NROOTS))
                {
                    uint di = loc[j] - pad;
                    
                    data[di] ^= cor;
                }
                else if (loc[j] < NN)
                {
                    uint pi = loc[j] - (NN - NROOTS);

                    parity[pi] ^= cor;
                }
            }
        }

        finish:

        if (erasPos != null)
        {
            for (int i = 0; i < count; i++)
                erasPos[i] = loc[i] - pad;
        }

        return count;
    }
}