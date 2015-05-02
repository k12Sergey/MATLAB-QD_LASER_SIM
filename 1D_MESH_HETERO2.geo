﻿// инициализация переменных
b1 = 0.2;
b2 = 0.3;
b3 = 0.2;


cl_1 = 0.01;
cl_2 = 0.001;


// ГРаница
Point(1) = {0, 0, 0, cl_1};
// гетеропереход
Point(2) = {b1, 0, 0, cl_2}; 
// гетеропереход
Point(3) = {b1+b2, 0, 0 ,cl_2};
// ГРаница
Point(4) = {b1+b2+b3, 0, 0, cl_1};

// линии 
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
 Line Loop(5) = {1,2,3};

// физические сущности
Physical Line (2001) = {1};
Physical Line (2002) = {2};
Physical Line (2003) = {3};  