# Решение системы трёхдиагональных матриц итерационными методами и методом прогонки
Релизация метода прогонки и итерационных методов: Якоби, Зейделя, Наискорейшего Спуска и Верхней релаксации для трёхдиагональной матрицы и решение ими системы линейных алгебраических уравнений:
$$(a_1+a_2+h^2g_1)y_1-a_2y_2=f_1h^2$$
$$\dots\dots\dots\dots\dots\dots$$
 $$-a_iy_{i-1}+(a_i+a_{i+1}+h^2g_i)y_i-a_{i+1}y_{i+1}=f_ih^2$$
 $$\quad\quad\quad \dots\dots\dots\dots\dots\dots $$
  $$(a_{n-1}+a_n+h^2g_{n-1})y_{n-1}-a_{n-1}y_{n-2}=f_{n-1}h^2$$
Здесь: $a_i=p(ih)$, $g_i=q(ih)$, $f_i=f(ih)$,
$f(x)=-(p(x)u'(x))'+q(x)u(x)$,
$h=\displaystyle\frac{1}{n}$

$u(x)=x^\alpha(1-x)^\beta;$

$p(x)=1+x^\gamma;$

$g(x)=x+1.$

Определение оптимального параметра $\omega$ метода верхней релаксации для заданной системы, а также вычисление необходимого числа итераций.

Данные для таблиц погрешности выводятся в виде кода TeX.

Программа также содержит класс трёхдиагональных матриц, методы создания таблиц TeX для любых двумерных массивов строк, CSV файлов, а также реализацию необходимых для работы алгебраических операций.
