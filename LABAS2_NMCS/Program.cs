
using System.Diagnostics;

namespace LABAS2_NMCS
{
    internal class Program
    {
        public class LUDecomposition
        {
            private readonly double[,] L;
            private readonly double[,] U;

            public LUDecomposition(double[,] A)
            {
                int n = A.GetLength(0);  
                L = new double[n, n];   
                U = new double[n, n];   

                // прямий хід метода Гаусса: знаходим U
                for (int i = 0; i < n; i++)
                {
                    L[i, i] = 1.0;  // елементи діагонали матриці L = 1
                    for (int j = i; j < n; j++)
                    {
                        double sum = 0.0;  // сумма для знахождения елемента U[i, j]
                        for (int k = 0; k < i; k++)
                        {
                            sum += L[i, k] * U[k, j];  // віднімаємо добуток елементів матриць L та U
                        }
                        U[i, j] = A[i, j] - sum;  // знаходим елемент матриці U
                    }

                    // зворотний хід методу Гауса: знаходимо L
                    for (int j = i + 1; j < n; j++)
                    {
                        double sum = 0.0;  // сумма для знахождення елемента L[j, i]
                        for (int k = 0; k < i; k++)
                        {
                            sum += L[j, k] * U[k, i];  // віднімаємо добуток елементів матриць L та U
                        }
                        L[j, i] = (A[j, i] - sum) / U[i, i];  // знаходимо елемент матриці L
                    }
                }
            }
            public double[] Solve(double[] b)
            {
                int n = b.Length;
                double[] y = new double[n];
                double[] x = new double[n];

                // На першому етапі знаходимо проміжний вектор Y, використовуючи пряму підстановку: Ly = b
                for (int i = 0; i < n; i++)
                {
                    y[i] = b[i];

                    for (int j = 0; j < i; j++)
                    {
                        y[i] = y[i] - L[i, j] * y[j];
                    }
                }

                // на другому етапі знаходимо безпосередньо вектор розв’язків Х, застосовуючи зворотну підстановку: Ux = y
                for (int i = n - 1; i >= 0; i--)
                {
                    x[i] = y[i];

                    for (int j = i + 1; j < n; j++)
                    {
                        x[i] = x[i] - U[i, j] * x[j];
                    }

                    x[i] = x[i] / U[i, i];
                }

                return x;
            }
        }

        static void Main(string[] args)
        {
            Console.WriteLine("Введіть 1 для того щоб розв'язати методом LU розкладу, або введіть 2, щоб розв'язати методом простої ітерації");
            int flag = Convert.ToInt32(Console.ReadLine());

            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            if (flag == 1)
            {
                 //double[,] A = { { -2, -3, 11 }, { 1, 12, -5 }, { -1, 1, -1 } };
                 //double[,] duplicate_A = { { -2, -3, 11 }, { 1, 12, -5 }, { -1, 1, -1 } };
                 //double[] b = { -15, 40, 35 };

                 double[,] A = { { -2, 1, 1 }, { 1, -2, 1 }, { -1, 3, -6 } }; 
                 double[,] duplicate_A = { { -2, 1, 1 }, { 1, -2, 1 }, { -1, 3, -6 } }; 
                  double[] b = { 15, 10, 12 };

                LUDecomposition lu = new LUDecomposition(A);
                double[] x = lu.Solve(b);

                Console.WriteLine("Solution: [{0}, {1}, {2}]", x[0], x[1], x[2]);
                Console.WriteLine("Check b[1]: {0}", duplicate_A[0, 0] * x[0] + duplicate_A[0, 1] * x[1] + duplicate_A[0, 2] * x[2]);
                Console.WriteLine("Check b[2]: {0}", duplicate_A[1, 0] * x[0] + duplicate_A[1, 1] * x[1] + duplicate_A[1, 2] * x[2]);
                Console.WriteLine("Check b[3]: {0}", duplicate_A[2, 0] * x[0] + duplicate_A[2, 1] * x[1] + duplicate_A[2, 2] * x[2]);


            }
            else if (flag == 2)
            {
                double[,] A = { { -2, 1, 1 }, { 1, -2, 1 }, { -1, 3, -6 } };
                double[] b = { 15, 10, 12 };

                //double[,] A = { { -2, -3, 11 }, { 1, 12, -5 }, { -1, 1, -1 } };
                //double[] b = { -15, 40, 35 };
                double[] vec = { -90, 15, 8 };
                double tolerance = 0.1;
                int maxIterations = 10000;

                double[] x = SolveBySimpleIteration(A, b, vec, tolerance, maxIterations);

                // Print solution
                for (int i = 0; i < x.Length; i++)
                {
                    Console.WriteLine("x[{0}] = {1}", i, x[i]);
                }
            }
            stopwatch.Stop();
            TimeSpan ts = stopwatch.Elapsed;

            Console.WriteLine("Час виконання програми: {0:hh\\:mm\\:ss\\:fff}", ts);

        }


        public static double[] SolveBySimpleIteration(double[,] A, double[] b, double[] veс, double tolerance, int maxIterations)
        {
            int n = b.Length;
            double[] x0 = new double[n];
            double[] x = new double[n];
            double error = double.MaxValue;
            int iteration = 0;
            
            for (int i = 0; i < n; i++)
            {
                x0[i] = veс[i];
            }


            // Виконуєм ітерації до збіжності або максимальної кількості ітерацій
            while (error > tolerance && iteration < maxIterations)
            {
                //  Обчислюємо  x_k+1 за допомогою x_k
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j)
                        {
                            sum += A[i, j] * x0[j];
                        }
                    }
                    x[i] = (b[i] - sum) / A[i, i];
                }

                // Помилка обчислень
                error = 0;
                for (int i = 0; i < n; i++)
                {
                    error += Math.Abs(x[i] - x0[i]);
                }

                // Оновлюємо х0
                Array.Copy(x, x0, n);

                iteration++;
            }

            if (iteration == maxIterations && error > tolerance)
            {
                throw new Exception("Метод не зійшовся");
            }

            return x;
        }
    }
}