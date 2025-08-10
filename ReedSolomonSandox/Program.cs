using ReedSolomonSandox;

byte[] bytes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

Console.WriteLine("Original:  " + string.Join(", ", bytes));

int res = new ReedSolomon63().Decode(
    bytes,
    36,
    null,
    null!,
    0,
    null!);

Console.WriteLine("Corrected: " + string.Join(", ", bytes));

Console.WriteLine(res);