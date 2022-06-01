using System;
using System.Collections.Generic;
using System.IO;

namespace ThomasAlgorithmSharp
{
    /// <summary>
    /// Класс для вывода данных в отчёт в виде кода TeX и создания CSV файлов
    /// </summary>
    class Reporter
    {
        public Reporter() { }
        public void Add(string nameTable, string[] names, string[][] data)
        {
            Console.WriteLine(PrintTableLaTeX(nameTable, names, data));
            //CreateCSV(nameTable, names, data);
        }
        /// <summary>
        /// Выводит в консоль данные в виде кода для таблицы TeX 
        /// </summary>
        /// <param name="nameTable"> Название таблицы </param>
        /// <param name="names"> Шапка таблицы </param>
        /// <param name="data"> Содержимое </param>
        /// <returns></returns>
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
            StreamWriter writer = new StreamWriter(nameFile + ".csv");
            for (int i = 0; i < names.Length - 1; i++)
            {
                writer.Write(names[i] + ";");
            }
            writer.WriteLine(names[names.Length - 1]);
            int rows = data[0].Length;
            int columns = data.Length;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns - 1; j++)
                {
                    writer.Write(data[j][i] + ";");
                }
                writer.WriteLine(data[columns - 1][i]);
            }
            writer.Close();
        }
    }
}
