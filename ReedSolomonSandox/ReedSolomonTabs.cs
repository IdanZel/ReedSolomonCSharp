namespace ReedSolomonSandox;

public class ReedSolomonTabs
{
    public const int SYM = 6;
    public const uint PRM = 1;
    public const uint MM = SYM;
    public const uint SIZE = 63; // maximum symbols in field
    public const uint NN = SIZE;
    public const uint A0 = SIZE;
    public const uint MODS = 0; // modulo table: 1/2 the symbol size squared, up to 4k
    public readonly uint Iprim; // 0 if uninitialized

    protected readonly byte[] AlphaTo = new byte[NN + 1];
    protected readonly byte[] IndexOf = new byte[NN + 1];

    public ReedSolomonTabs()
    {
        // Do init if not already done.  We check one value which is initialized to 0; this is
        // safe, 'cause the value will not be set 'til the initializing thread has completely
        // initialized the structure.  Worst case scenario: multiple threads will initialize
        // identically.  No mutex necessary.
        if (Iprim > 0)
            return;

        // Generate Galois field lookup tables
        IndexOf[0] = (byte)A0; // log(zero) = -inf
        AlphaTo[A0] = 0; // alpha**-inf = 0
        uint sr = GaloisFieldPoly.Op(0);
        for (int i = 0; i < NN; i++)
        {
            IndexOf[sr] = (byte)i;
            AlphaTo[i] = (byte)sr;
            sr = GaloisFieldPoly.Op(sr);
        }

        // If it's not primitive, raise exception or abort
        if (sr != AlphaTo[0])
        {
            throw new Exception("reed-solomon: Galois field polynomial not primitive");
        }
        
        // Find prim-th root of 1, index form, used in decoding.
        uint iptmp = 1;
        while (iptmp % PRM != 0)
            iptmp += NN;
        Iprim = iptmp / PRM;
    }

    protected static byte Modnn(uint x)
    {
        while (x >= NN + MODS)
        {
            x -= NN;
            x = (x >> (int)MM) + (x & NN);
        }

        return (byte)x;
    }
}