syms h n x
limit((cos(x+h)-cos(x))/h, h, 0)
limit((1 + x/n)^n, n, inf)

%% Linksseitig
syms x
limit(x/abs(x), x, 0, "left")

%% Rechtsseitig
syms x
limit(x/abs(x), x, 0, "right")

%% Beidseitig
syms x
limit (x/abs(x), x, 0)

%% Function test
syms f x
f = x/abs(x)

[left, right, both] = compute_limits(f, 0);

%% Ü1 Stetigkeit (Keine funktionierende Lösung)
syms a x f

f = piecewise(x<2, 8*a+16*x, x>= 2, a^2*(x+2));
x0 = 2;

ermittleStetigkeitDual(f, x0, a)
%Funktioniert nicht, da a seperat betrachtet werden muss?

%% Ü2 Stetigkeit (Keine funktionierende Lösung)
syms a x f
f = piecewise(x<=2, 2*x-6, x>2, a^2-a*x^2+2);
x0 = 2;

ermittleStetigkeitDual(f, x0, a)

%% Ü3 Stetigkeit (funktionierende Lösung)
syms x f1 f2
f1 = x^2+1;
f2 = -x*(x-3);

ermittleStetigkeit(f1, f2, 1, true)

%% Test Differenzierbarkeit
syms x f;

f = piecewise(x<=1, x^2+1, x>1, -x*(x-3));
x0=1;
ermittleDifferenzierbarkeit(f, x0)

%% Test Master-Theorem Stefan
syms n;

a=8;
b=2;
f = n^2+7*n+5;
logb_a = log(a) / log(b);

if simplify(f / n^(logb_a - 1)) == 0
    disp("1. Methode")
elseif simplify(f /n^logb_a) == 1
    disp("2. Methode")
elseif simplify (f /n^(logb_a+1)) == 1
    disp("3. Methode")
else
    disp("Keine Methode anwendbar")
end

%% Test Master-Theorem
syms f n;
a = 8;
b = 2;
f = n^2+7*n+5;
masterTheorem(a, b, f)
%% Function Limits (Working)
function [left, right, both] = compute_limits(f, x0)
    syms x
    left = limit(f, x, x0, "left");
    right = limit(f, x, x0, "right");
    both = limit(f, x, x0);

    disp("Links:")
    disp(left)
    disp("Rechts:")
    disp(right)
    disp("Beide:")
    disp(both)
end

%% function Stetigkeit
function [istStetig] = ermittleStetigkeit(f, x0)
    syms x;
    istStetig = limit(f, x, x0) == subs(f, x, x0);
    if istStetig
        disp("FUnktion ist an x0 = "+string(x0) +" stetig.")
    else
        disp("Funktion ist an x0 = "+string(x0) +" diskret.")
    end
end

%% function Stetigkeit mit Zweitvariable a (Not working)
function [istStetig, a_vals] = ermittleStetigkeitDual(f, x0, zweitVar)
    syms x;
    istStetig = limit(f, x, x0, "left") == limit(f, x, x0, "right");
    if istStetig
        disp("Funktion f ist stetig an x0 = "+string(x0))
    else
        disp("Funktion f ist diskret an x0 = "+string(x0))
    end
  
    a_vals = solve(limit(f, x, x0, "left")==limit(f, x, x0, "right"), zweitVar);
    disp("Die Werte von "+string(zweitVar)+" bei denen die Funktion an x0 = "+string(x0)+" stetig sind:")
    disp(a_vals)
end

%% function Differenzierbarkeit
function istDiffbar = ermittleDifferenzierbarkeit(f, x0)
    syms x;
    f_diff = diff(f, x);
    
    left = limit(f_diff, x, x0, "left");
    right = limit(f_diff, x, x0, "right");
    istDiffbar = left == right;
    if istDiffbar
        disp("Funktion f ist an x0 = "+string(x0)+" differenzierbar!")
    else
        disp("Funktion f ist an x0 = "+string(x0)+" nicht differenzierbar!")
    end
end

%% Master-Theorem
function masterTheorem(a, b, f)
    syms n;
    logb_a = log(a)/log(b);

    if simplify(f/n^(logb_a-1)) == 0
        disp("Gemäß dem Master-Theorem: T(n) = Theta(n^log_b(a))");
    elseif simplify(f / n^(logb_a))==1
        disp("Gemäß dem Master-Theorem: T(n) = Theta(n^log_b(a) * log(n))");
    elseif simplify(f/n^(logb_a+1))==0
        disp("Gemäß dem Master-Theorem: T(n) = Theta(f(n))");
    else
        disp("Master-Theorem nicht anwendbar!")
    end
end
    