using System;

namespace ThomasAlgorithmSharp
{
    class Program
    {
        static double alpha =2.0;
        static double beta =1.0;
        static double gamma =1.0;
        static int n =10;
        static double h = 1.0/n;
        static int N = n+1;
        static double epsilon = Math.Pow(h, 3);
        #region functions
        static double Q(double x)
        {
            return x + 1;
        }
        static double P(double x)
        {
            return 1 + Math.Pow(x, gamma);
        }
        static double U(double x)
        {
            return Math.Pow(x, alpha) * Math.Pow(1.0 - x, beta);
        }
        static double a_i(int i)
        {
            return P(i * h);
        }
        static double q_i(int i)
        {
            return Q(i * h);
        }
        static double F(int i) //qu-p'u'-pu''
        {
            double x = i * h;
            return Q(x) * U(x) -
                gamma * Math.Pow(x, gamma - 1) * //p'
                (alpha * Math.Pow(x, alpha - 1) * Math.Pow(1 - x, beta) - beta * Math.Pow(1 - x, beta - 1) * Math.Pow(x, alpha)) - //u'
                P(x) * ( //p * u''
                alpha * (alpha - 1) * Math.Pow(x, alpha - 2) * Math.Pow(1 - x, beta) -
                alpha * beta * Math.Pow(1 - x, beta - 1) * Math.Pow(x, alpha - 1) +
                beta * (beta - 1) * Math.Pow(1 - x, beta - 2) * Math.Pow(x, alpha) -
                alpha * beta * Math.Pow(x, alpha - 1) * Math.Pow(1 - x, beta - 1)
                );
        }
        #endregion

        #region r
        static double[] R(DiagonalMatrix a, double[] y, double[] b)
        {
            double[] r = new double[y.Length];
            for (int i = 0; i < r.Length; i++)
            {
                for (int j = 0; j < r.Length; j++)
                {
                    r[i] += a[i, j] * y[j];
                }
                r[i] -= b[i];
            }
            return r;
        }
        static double MaxR(DiagonalMatrix a, double[] y, double[] b)
        {
            double[] r = R(a, y, b);
            double max = 0;
            foreach (double i in r)
            {
                if (max < Math.Abs(i)) max = Math.Abs(i);
            }
            return max;
        }
        #endregion

        #region initiate
        static DiagonalMatrix InitiateMatrix()
        {
            DiagonalMatrix matrix = new DiagonalMatrix(N);
            matrix.SetDiag(0, 1);
            for (int i = 1; i < n-1; i++)
            {
                matrix.SetHighDiag(i, -P((i + 1) * h));
            }
            for (int i = 2; i < n; i++)
            {
                matrix.SetLowDiag(i, -P(i * h));
            }
            for (int i = 1; i < n; i++)
            {
                matrix.SetDiag(i, P(i * h) + P((i + 1) * h) + h * h * Q(i * h));
            }
            matrix.SetDiag(n, 1);
            return matrix;
        }
        static double[] InitiateB()
        {
            double[] result = new double[N];
            result[0] = 0;
            for (int i = 1; i < n; i++)
            {
                if (i == 0)
                {
                    result[i] = 0;
                }
                else if (i == n)
                {
                    result[i] = 0;
                }
                else
                {
                    result[i] = F(i) * h * h;
                }
            }
            result[n] = 0;
            return result;
        }

        #endregion

        #region ThomasAlgorithm
        static void DirectCourse(DiagonalMatrix matrix, double[] b)
        {
            double a;
            for (int i = 2; i<n; i++)
            {
                a = matrix[i, i - 1] / matrix[i - 1, i - 1];
                matrix[i, i - 1] -= matrix[i - 1, i - 1] * a;
                matrix[i, i] -= matrix[i - 1, i] * a;
                b[i] -= b[i - 1] * a;
            }
            matrix.Print(b);
        }
        static double[] ReverseCourse(DiagonalMatrix matrix, double[] b)
        {
            double[] y = new double[b.Length];
            y[n] = 0;
            for (int i = n - 1; i > 0; i--)
            {
                y[i] = (b[i]  - y[i + 1] * matrix[i, i + 1])/ matrix[i, i];
            }
            return y;
        }
        static void Thomas(DiagonalMatrix A)
        {
            DiagonalMatrix matrix = A.Copy();
            double[] B = InitiateB();
            matrix.Print(B);
            DirectCourse(matrix, B);
            double[] y = ReverseCourse(matrix, B);
            Console.WriteLine("Метод Прогонки");
            PrintTable(y);
        }
        #endregion
        static void Seidel(DiagonalMatrix A)
        {
            DiagonalMatrix a = A.Copy();
            double[] b = InitiateB();
            double[] x = new double[N];
            int k = 0;
            do
            {
                for (int i = 0; i < x.Length; i++)
                {
                    x[i] = b[i] / a[i, i];
                    for (int j = 0; j < x.Length; j++)
                    {
                        if (i != j)
                        {
                            x[i] -= a[i, j] / a[i, i] * x[j];
                        }
                    }
                }
                k++;
            }
            while (MaxR(a, x, b) > epsilon);
            Console.WriteLine($"Метод Зейделя {k} итераций");
            PrintTable(x);
        }
        static void Jacobi(DiagonalMatrix A)
        {
            DiagonalMatrix a = A.Copy();
            double[] b = InitiateB();
            double[] x = new double[N];
            double[] cx = new double[N];
            int k = 0;
            do
            {
                for (int i = 0; i < x.Length; i++)
                {
                    cx[i] = b[i] / a[i, i];
                    for (int j = 0; j < x.Length; j++)
                    {
                        if (i != j)
                        {
                            cx[i] -= a[i, j] / a[i, i] * x[j];
                        }
                    }
                }
                for (int i = 0; i < x.Length; i++)
                {
                    x[i] = cx[i];
                }
                k++;
            }
            while (MaxR(a, x, b) > epsilon);
            Console.WriteLine($"Метод Якоби {k} итераций");
            PrintTable(x);
        }
        static int Relax(DiagonalMatrix A, double omega)
        {
            DiagonalMatrix a = A.Copy();
            double[] b = InitiateB();
            double[] x = new double[N];
            double[] cx = new double[N];
            int k = 0;
            do
            {
                for (int i =0; i<x.Length; i++)
                {
                    cx[i] = b[i] / a[i, i];
                    //cx[i] = (1 - omega) * x[i];
                    for (int j = 0; j<i; j++)
                    {
                        cx[i] -= a[i, j] / a[i, i] * cx[j];
                    }
                    for (int j =i+1; j < x.Length; j++)
                    {
                        cx[i] -= a[i, j] / a[i, i] * x[j];
                    }
                    cx[i] *= omega;
                    cx[i] += (1 - omega) * x[i];
                }
                for (int i = 0; i < x.Length; i++)
                {
                    x[i] = cx[i];
                }
                k++;
            }
            while (MaxR(a, x, b) > epsilon);
            for (int i = 0; i < x.Length; i++)
            {
                Console.WriteLine(x[i]);
            }
            return k;
        }
        static void PrintTable(double[] x)
        {
            Console.WriteLine("{0,20}{1,25}{2,30}{3,30}", "ih", "yi", "u(ih)", "|yi-u(ih)|");
            for (int i =0; i<x.Length; i++)
            {
                double u = U(i * h);
                Console.WriteLine("{0,20}{1,25}{2,30}{3,30}", i * h, x[i], u, Math.Abs(x[i] - u));
            }

        }
        static void Main(string[] args)
        {
            DiagonalMatrix A = InitiateMatrix();
            //PrintTable();
            Thomas(A);
            Seidel(A);
            //A.Print();
            Jacobi(A);
        }
        

    }
    class DiagonalMatrix
    {
        double[] diag;
        double[] high_diag;
        double[] low_diag;
        int size;
        public DiagonalMatrix(int n)
        {
            size = n;
            diag = new double[n];
            high_diag = new double[n - 1];
            low_diag = new double[n - 1];
        }
        public DiagonalMatrix Copy()
        {
            DiagonalMatrix clone = new DiagonalMatrix(size);
            for (int i = 0; i < high_diag.Length; i++)
            {
                clone.SetHighDiag(i, high_diag[i]);
            }
            for (int i = 0; i < diag.Length; i++)
            {
                clone.SetDiag(i, diag[i]);
            }
            for (int i = 0; i < low_diag.Length; i++)
            {
                clone.SetLowDiag(i+1, low_diag[i]);
            }
            return clone;
        }
        public int GetLength()
        {
            return size;
        }
        public double this[int row, int col]
        {
            get
            {
                if (row == col)
                {
                    return diag[row];
                }
                else if (col - row == 1)
                {
                    return high_diag[row];
                }
                else if (row - col == 1)
                {
                    return low_diag[col];
                }
                else return 0;
            }
            set
            {
                if (row == col)
                {
                    diag[row] = value;
                }
                else if (col - row == 1)
                {
                    high_diag[row] = value;
                }
                else if (row - col == 1)
                {
                    low_diag[col] = value;
                }
                else
                {
                    Console.WriteLine("Попытка записать не диагональный элемент");
                }
            }
        }
        public void SetHighDiag(int row, double value)
        {
            high_diag[row] = value;
        }
        public void SetDiag(int row, double value)
        {
            diag[row] = value;
        }
        public void SetLowDiag(int row, double value)
        {
            low_diag[row-1] = value;
        }
        public void Print()
        {
            Console.WriteLine("Матрица коэффициентов:");
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    Console.Write(this[i,j].ToString()+' ');
                }
                Console.WriteLine();
            }
        }
        public void Print(double[] b)
        {
            Console.WriteLine("Матрица коэффициентов | значение:");
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    Console.Write(this[i, j].ToString() + ' ');
                }
                Console.WriteLine(" | "+b[i]);
            }
        }
        static public void TestMatrix(DiagonalMatrix matrix)
        {
            for (int i = 0; i < matrix.GetLength(); i++)
            {
                for (int j = 0; j < matrix.GetLength(); j++)
                {
                    if (i == j)
                    {
                        matrix[i, j] = 2;
                    }
                    else if (j - i == 1)
                    {
                        matrix[i, j] = 3;
                    }
                    else if (i - j == 1)
                    {
                        matrix[i, j] = 1;
                    }
                }
            }
            matrix.Print();
        }
    }
}
