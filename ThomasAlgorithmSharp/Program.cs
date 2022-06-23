using System;
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
        static double epsilon = Math.Pow(h, 3); //задание величины эпсилон-окрестности
        static string[] headers = { "ih", "yi", "u(ih)", "|yi-u(ih)|" }; //шапка таблиц погрешности
        static double bestOmega; //лучший параметр метода релаксации для заданной системы
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
        /// <summary>
        /// Вычисление компонент вектора невязки r = Ay-b
        /// </summary>
        /// <param name="A"> Матрица </param>
        /// <param name="y"> Вектор значений(в виде массива) </param>
        /// <param name="b"> Вектор правой части(в виде массива) </param>
        /// <returns> Массив компонент вектора невязки </returns>
        static double[] R(TridiagonalMatrix A, double[] y, double[] b)
        {
            double[] r = new double[y.Length];
            for (int i = 0; i < r.Length; i++)
            {
                for (int j = 0; j < r.Length; j++)
                {
                    r[i] += A[i, j] * y[j];
                }
                r[i] -= b[i];
            }
            return r;
        }
        /// <summary>
        /// Вычисляет вектор невязки и возвращает его максимальную компоненту
        /// </summary>
        /// <param name="A"> Матрица </param>
        /// <param name="y"> Вектор значений(в виде массива) </param>
        /// <param name="b"> Вектор правой части(в виде массива) </param>
        /// <returns> Возвращает максимальную по модулю компонету вектора невязки </returns>
        static double MaxR(TridiagonalMatrix A, double[] y, double[] b)
        {
            double[] r = R(A, y, b);
            double max = 0;
            foreach (double i in r)
            {
                if (max < Math.Abs(i)) max = Math.Abs(i);
            }
            return max;
        }
        #endregion

        #region algem
        static double[] MultiplyingMatrixVector(TridiagonalMatrix A, double[] x)
        {
            double[] res = new double[x.Length];
            for (int i = 0; i<x.Length; i++)
            {
                for (int j = 0; j<x.Length; j++)
                {
                    res[i] += A[i, j] * x[j];
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
        /// <summary>
        /// Создание матрицы системы(система описана в файле README.md)
        /// </summary>
        /// <returns> Матрица левой части системы </returns>
        static TridiagonalMatrix InitiateMatrix()
        {
            TridiagonalMatrix matrix = new TridiagonalMatrix(N);
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
        /// <summary>
        /// Создание вектора правой части системы(система описана в файле README.md)
        /// </summary>
        /// <returns> Массив компонент правой части </returns>
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

        /// <summary>
        /// Нахождение погрешности в точке
        /// </summary>
        /// <param name="y"> Значение </param>
        /// <param name="i"> Номер точки </param>
        /// <returns> абсолютная погрешность в точке </returns>
        static double Inaccuracy(double y, int i)
        {
            return Math.Abs(y - U(i * h));
        }

        #region ThomasAlgorithm
        /// <summary>
        /// Прямой ход метода Прогонки
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="b"></param>
        static void DirectCourse(TridiagonalMatrix matrix, double[] b)
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
        /// <summary>
        /// Обратный ход метода Прогонки
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        static double[] ReverseCourse(TridiagonalMatrix matrix, double[] b)
        {
            double[] y = new double[b.Length];
            y[n] = 0;
            for (int i = n - 1; i > 0; i--)
            {
                y[i] = (b[i]  - y[i + 1] * matrix[i, i + 1])/ matrix[i, i];
            }
            return y;
        }
        /// <summary>
        /// Метод Прогонки
        /// </summary>
        /// <param name="A"></param>
        static void Thomas(TridiagonalMatrix A)
        {
            TridiagonalMatrix matrix = (TridiagonalMatrix)A.Clone();
            double[] B = InitiateB();
            DirectCourse(matrix, B);
            double[] y = ReverseCourse(matrix, B);
            string header = "Метод Прогонки";
            reporter.Add(header, headers, DataToTable(y));
        }
        #endregion
        /// <summary>
        /// Метод Зейделя
        /// </summary>
        /// <param name="A"></param>
        static void Seidel(TridiagonalMatrix A)
        {
            TridiagonalMatrix a = (TridiagonalMatrix)A.Clone();
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

        /// <summary>
        /// Метод Якоби
        /// </summary>
        /// <param name="A"></param>
        static void Jacobi(TridiagonalMatrix A)
        {
            TridiagonalMatrix a = (TridiagonalMatrix)A.Clone();
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
        /// <summary>
        /// Поиск оптимального параметра омега для метода верхней релаксации
        /// </summary>
        /// <param name="A"></param>
        /// <param name="numberParts"> Число делений отрезка [1;2]</param>
        static void RelaxFromOmega(TridiagonalMatrix A, int numberParts)
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
        /// <summary>
        /// Метод релаксации с заданным параметром omega
        /// </summary>
        /// <param name="A"></param>
        /// <param name="omega"></param>
        /// <returns> Число итераций, необходимых для решения при заданном параметре </returns>
        static int Relax(TridiagonalMatrix A, double omega)
        {
            TridiagonalMatrix a = (TridiagonalMatrix)A.Clone();
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

       /// <summary>
       /// Метод Верхней Релаксации
       /// </summary>
       /// <param name="A"></param>
        static void Relax(TridiagonalMatrix A)
        {
            TridiagonalMatrix a = (TridiagonalMatrix)A.Clone();
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
        #region tables
        /// <summary>
        /// Формирование данных для таблицы погрешности
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
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
        /// <summary>
        /// Формирование таблицы из двух столбцов: точка, значение
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
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
        #endregion tables
        /// <summary>
        /// Метод Наискорейшего Спуска
        /// </summary>
        /// <param name="A"></param>
        static void FastestDescent(TridiagonalMatrix A)
        {
            TridiagonalMatrix a = (TridiagonalMatrix)A.Clone();
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
            TridiagonalMatrix A = InitiateMatrix();
            Thomas(A);
            Jacobi(A);
            Seidel(A);
            FastestDescent(A);
            RelaxFromOmega(A, 20);
            Relax(A);
        }
    }
}
