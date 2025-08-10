namespace ReedSolomonSandox;

public static class GaloisFieldPoly
{
    private const int SYM = 6;
    private const uint POLY = 0x43;

    public static uint Op(uint sr)
    {
        if (sr == 0)
            sr = 1;
        else
        {
            sr <<= 1;
            if ((sr & (1 << SYM)) > 0)
                sr ^= POLY;
            sr &= ((1 << SYM) - 1);
        }

        return sr;
    }
}