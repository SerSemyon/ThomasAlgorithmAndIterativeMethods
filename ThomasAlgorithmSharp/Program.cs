using System;
using System.IO;
using System.Collections.Generic;
using System.Collections;
namespace ThomasAlgorithmSharp
{
    class Program
    {
        static double alpha = 1.0;
        static double beta = 1.0;
        static double gamma = 2.0;
        static int n = 10;
        static double h = 1.0 / n;
        static int N = n + 1;
        static double epsilon = Math.Pow(h, 3);
        static string[] headers = { "ih", "yi", "u(ih)", "|yi-u(ih)|" };
        static double bestOmega;
        static Reporter reporter = new Reporter();
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
                (alpha * Math.Pow(x, alpha - 1) * Math.Pow(1 - x, beta) - 
                beta * Math.Pow(1 - x, beta - 1) * Math.Pow(x, alpha)) - //u'
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

        #region algem
        static double[] MultiplyingMatrixVector(DiagonalMatrix a, double[] x)
        {
            double[] res = new double[x.Length];
            for (int i = 0; i<x.Length; i++)
            {
                for (int j = 0; j<x.Length; j++)
                {
                    res[i] += a[i, j] * x[j];
                }
            }
            return res;
        }
        static double ScalarProduct(double[] x, double[] y)
        {
            double res = 0.0;
            for (int i = 0; i<x.Length; i++)
            {
                res += x[i] * y[i];
            }
            return res;
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
            DirectCourse(matrix, B);
            double[] y = ReverseCourse(matrix, B);
            string header = "Метод Прогонки";
            reporter.Add(header, headers, DataToTable(y));
        }
        #endregion
        static void Seidel(DiagonalMatrix A)
        {
            DiagonalMatrix a = A.Copy();
            double[] b = InitiateB();
            double[] x = new double[N];
            List<double> inaccuracy = new List<double>();
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
                double cInaccuracy = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    if (cInaccuracy < Inaccuracy(x[i], i)) cInaccuracy = Inaccuracy(x[i], i);
                }
                k++;
                inaccuracy.Add(cInaccuracy);
            }
            while (MaxR(a, x, b) > epsilon);
            string header = $"Метод Зейделя. Число итераций {k}";
            string[] Headers = { "Число итераций", "Погрешность" };
            double[] vecInaccuracy = new double[inaccuracy.Count];
            double[] K = new double[inaccuracy.Count];
            for (int i = 0; i < vecInaccuracy.Length; i++)
            {
                K[i] = i + 1;
                vecInaccuracy[i] = inaccuracy[i];
            }
            reporter.Add(header, headers, DataToTable(x));
            header = "Метод Зейделя. График погрешности";
            Reporter.CreateCSV(header, Headers, Table2Columns(K, vecInaccuracy));
        }

        static void Jacobi(DiagonalMatrix A)
        {
            DiagonalMatrix a = A.Copy();
            double[] b = InitiateB();
            double[] x = new double[N];
            double[] cx = new double[N];
            List<double> inaccuracy = new List<double>();
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
                double cInaccuracy = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    if (cInaccuracy < Inaccuracy(cx[i],i)) cInaccuracy = Inaccuracy(cx[i], i);
                    x[i] = cx[i];
                }
                k++;
                inaccuracy.Add(cInaccuracy);
            }
            while (MaxR(a, x, b) > epsilon);
            string header = $"Метод Якоби. Число итераций {k}";
            string[] Headers = { "Число итераций", "Погрешность" };
            double[] vecInaccuracy = new double[inaccuracy.Count];
            double[] K = new double[inaccuracy.Count];
            for (int i = 0; i<vecInaccuracy.Length; i++)
            {
                K[i] = i + 1;
                vecInaccuracy[i] = inaccuracy[i];
            }
            reporter.Add(header, headers, DataToTable(x));
            header = "Метод Якоби. График погрешности";
            Reporter.CreateCSV(header, Headers, Table2Columns(K, vecInaccuracy));
        }
        #region Relax
        static void RelaxFromOmega(DiagonalMatrix A, int numberParts)
        {
            double omegaH = 1.0 / numberParts;
            double[] Omega = new double[numberParts];
            double[] inaccuracy = new double[numberParts];
            Omega[0] = 1;
            double omega = 1;
            for (int i = 0; i < Omega.Length; i++)
            {
                Omega[i] = omega;
                inaccuracy[i] = Relax(A, omega);
                omega += omegaH;
            }
            string header = "Омега";
            string[] Headers = { @"\omega", "Число итераций" };
            string[][] data = Table2Columns(Omega, inaccuracy);
            reporter.Add(header, Headers, data);
            double optimOmega = Omega[0];
            double min = inaccuracy[0];
            for (int i = 0; i<Omega.Length; i++)
            {
                if (inaccuracy[i] < min)
                {
                    min = inaccuracy[i];
                    optimOmega = Omega[i];
                }
            }
            bestOmega = Math.Round(optimOmega, 2);
            Console.WriteLine("Оптимальное значение омега - " + bestOmega);
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
            return k;
        }

        static void Relax(DiagonalMatrix A)
        {
            DiagonalMatrix a = A.Copy();
            double[] b = InitiateB();
            double[] x = new double[N];
            double[] cx = new double[N];
            List<double> inaccuracy = new List<double>();
            int k = 0;
            do
            {
                for (int i = 0; i < x.Length; i++)
                {
                    cx[i] = b[i] / a[i, i];
                    for (int j = 0; j < i; j++)
                    {
                        cx[i] -= a[i, j] / a[i, i] * cx[j];
                    }
                    for (int j = i + 1; j < x.Length; j++)
                    {
                        cx[i] -= a[i, j] / a[i, i] * x[j];
                    }
                    cx[i] *= bestOmega;
                    cx[i] += (1 - bestOmega) * x[i];
                }
                double cInaccuracy = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    if (cInaccuracy < Inaccuracy(cx[i], i)) cInaccuracy = Inaccuracy(cx[i], i);
                    x[i] = cx[i];
                }
                k++;
                inaccuracy.Add(cInaccuracy);
            }
            while (MaxR(a, x, b) > epsilon);
            string header = $"Метод Релаксации. Число итераций {k}";
            string[] Headers = { "Число итераций", "Погрешность" };
            double[] vecInaccuracy = new double[inaccuracy.Count];
            double[] K = new double[inaccuracy.Count];
            for (int i = 0; i < vecInaccuracy.Length; i++)
            {
                K[i] = i + 1;
                vecInaccuracy[i] = inaccuracy[i];
            }
            reporter.Add(header, headers, DataToTable(x));
            header = "Метод Релаксации. График погрешности";
            Reporter.CreateCSV(header, Headers, Table2Columns(K, vecInaccuracy));
        }
        #endregion relax
        static void PrintTable(double[] x)
        {
            Console.WriteLine("{0,20}{1,25}{2,30}{3,30}", "ih", "yi", "u(ih)", "|yi-u(ih)|");
            for (int i =0; i<x.Length; i++)
            {
                double u = U(i * h);
                Console.WriteLine("{0,20}{1,25}{2,30}{3,30}", i * h, x[i], u, Math.Abs(x[i] - u));
            }

        }
        static string[][] DataToTable(double[] x)
        {
            string[][] data = new string[4][];
            for (int i = 0; i<4; i++)
            {
                data[i] = new string[x.Length];
            }
            for (int i =0; i<x.Length; i++)
            {
                data[0][i] = Math.Round(h * i,2).ToString();
                data[1][i] = Math.Round(x[i],10).ToString();
                double u = U(i * h);
                data[2][i] = Math.Round(u,10).ToString();
                data[3][i] = Math.Round(Math.Abs(x[i] - u),10).ToString();
            }
            return data;
        }
        static string[][] Table2Columns(double[] x, double[] y)
        {
            string[][] data = new string[2][];
            data[0] = new string[x.Length];
            data[1] = new string[y.Length];
            for (int i =0; i < x.Length; i++)
            {
                data[0][i] = Math.Round(x[i],3).ToString();
                data[1][i] = Math.Round(y[i],10).ToString();
            }
            return data;
        }
        static double Inaccuracy(double x, int i)
        {
            return Math.Abs(x - U(i * h));
        }
        static void FastestDescent(DiagonalMatrix A)
        {
            DiagonalMatrix a = A.Copy();
            double[] b = InitiateB();
            double[] x = new double[N];
            List<double> inaccuracy = new List<double>();
            int k = 0;
            double r;
            double[] rk;
            do
            {
                rk = R(a, x, b);
                r = ScalarProduct(rk, rk) / ScalarProduct(MultiplyingMatrixVector(a, rk), rk);
                for (int i = 0; i <x.Length; i++)
                {
                    x[i] -= r * rk[i];
                }
                double cInaccuracy = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    if (cInaccuracy < Inaccuracy(x[i], i)) cInaccuracy = Inaccuracy(x[i], i);
                }
                k++;
                inaccuracy.Add(cInaccuracy);
            }
            while (MaxR(a, x, b) > epsilon);
            string header = $"Метод наискорейшего спуска. Число итераций {k}";
            string[] Headers = { "Число итераций", "Погрешность" };
            double[] vecInaccuracy = new double[inaccuracy.Count];
            double[] K = new double[inaccuracy.Count];
            for (int i = 0; i < vecInaccuracy.Length; i++)
            {
                K[i] = i + 1;
                vecInaccuracy[i] = inaccuracy[i];
            }
            reporter.Add(header, headers, DataToTable(x));
            header = "Метод наискорейшего спуска. График погрешности";
            Reporter.CreateCSV(header, Headers, Table2Columns(K, vecInaccuracy));
        }
        static void Main(string[] args)
        {
            DiagonalMatrix A = InitiateMatrix();
            Thomas(A);
            Jacobi(A);
            Seidel(A);
            FastestDescent(A);
            RelaxFromOmega(A, 20);
            Relax(A);
        }
        static void AllMethod(int number)
        {
            n = number;
            N = n + 1;
        }
    }
    class Reporter
    {
        public void Add(string nameTable, string[] names, string[][] data)
        {
            Console.WriteLine(PrintTableLaTeX(nameTable, names, data));
            //CreateCSV(nameTable, names, data);
        }
        static string PrintTableLaTeX(string nameTable, string[] names, string[][] data)
        {
            //настройки таблицы
            string result = @"\begin{table}\caption{" + nameTable + @"}\begin {tabular}{|";
            for (int i = 0; i < names.Length; i++)
            {
                result += @"p{ 3cm}|";
            }
            result += @"} \hline $" + names[0] + "$ \n";
            for (int i = 1; i < names.Length; i++) //шапка таблицы
            {
                result += " & $" + names[i] + "$";
            }
            result += @"\\ \hline ";
            int rows = data[0].Length;
            int columns = data.Length;
            for (int i = 0; i < rows; i++)
            {
                if (i != 0) result += @"\\ \hline";
                result += "\n" + $"{data[0][i]}";
                for (int j = 1; j < columns; j++)
                {
                    result += $" & {data[j][i]}";
                }
            }
            result += @"\\ \hline\end{tabular} \end{table}" + "\n";
            return result;
        }
        public static void CreateCSV(string nameFile, string[] names, string[][] data)
        {
            StreamWriter writer = new StreamWriter("Данные для графиков\\"+nameFile + ".csv");
            for (int i = 0; i<names.Length-1; i++)
            {
                writer.Write(names[i]+";");
            }
            writer.WriteLine(names[names.Length - 1]);
            int rows = data[0].Length;
            int columns = data.Length;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j<columns-1; j++)
                {
                    writer.Write(data[j][i] + ";");
                }
                writer.WriteLine(data[columns-1][i]);
            }
            writer.Close();
        }
        public Reporter() {}
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
