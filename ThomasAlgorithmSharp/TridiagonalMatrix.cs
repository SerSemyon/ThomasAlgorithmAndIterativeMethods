using System;

namespace ThomasAlgorithmSharp
{
    /// <summary>
    /// Класс трёхдиагональных матриц
    /// </summary>
    class TridiagonalMatrix
    {
        double[] diag;
        double[] high_diag;
        double[] low_diag;
        int size;
        public int Size
        {
            get
            {
                return size;
            }
        }
        public TridiagonalMatrix(int n)
        {
            size = n;
            diag = new double[n];
            high_diag = new double[n - 1];
            low_diag = new double[n - 1];
        }
        public TridiagonalMatrix Copy()
        {
            TridiagonalMatrix clone = new TridiagonalMatrix(size);
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
                clone.SetLowDiag(i + 1, low_diag[i]);
            }
            return clone;
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
        /// <summary>
        /// Задание элемента верхней диагонали
        /// </summary>
        /// <param name="row"> Индекс строки </param>
        /// <param name="value"> Записываемое значение </param>
        public void SetHighDiag(int row, double value)
        {
            high_diag[row] = value;
        }
        /// <summary>
        /// Задание элемента центральной диагонали
        /// </summary>
        /// <param name="row"> Индекс строки </param>
        /// <param name="value"> Записываемое значение </param>
        public void SetDiag(int row, double value)
        {
            diag[row] = value;
        }
        /// <summary>
        /// Задание элемента нижней диагонали
        /// </summary>
        /// <param name="row"> Индекс строки </param>
        /// <param name="value"> Записываемое значение </param>
        public void SetLowDiag(int row, double value)
        {
            low_diag[row - 1] = value;
        }
        public void Print()
        {
            Console.WriteLine("Матрица коэффициентов:");
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    Console.Write(this[i, j].ToString() + ' ');
                }
                Console.WriteLine();
            }
        }
        /// <summary>
        /// Выводит матрицу в виде системы уравнений с вектором значений b
        /// </summary>
        /// <param name="b"> Вектор значений </param>
        public void Print(double[] b)
        {
            Console.WriteLine("Матрица коэффициентов | значение:");
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    Console.Write(this[i, j].ToString() + ' ');
                }
                Console.WriteLine(" | " + b[i]);
            }
        }
    }
}
