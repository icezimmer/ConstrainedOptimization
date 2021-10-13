function [x_min, f_min, elapsedTime, numSteps]= FrankWolfeEVAL(Q, q, P, xStart, eps, eps_ls, options1, options2)
%{
Funzione che calcola il minimo
INPUT:
Q matrice nxn semidefinita positiva
q vettore
P matrice kxn: la riga k-esima di P indica l'insieme Ik (P(k,j) = 1 se j sta in Ik,
0 altrimenti)
xStart: vettore x di partenza appartenente al dominio
%}

if isequal(options1,'NM')
    disp('LineSearch Newton Method')
elseif isequal(options1,'LBM')
    disp('LineSearch Linear Bisection Method')
elseif isequal(options1,'QBM')
    disp('LineSearch Quadratic Bisection Method')
else
    disp('LineSearch Trivial Method')
end

if isequal(options2,'yes')
    disp('Using MOMENTUM')
end

%forzo q e x a essere vettori colonna
q = q(:);
x = xStart(:);

%fuzione f
f = @(x) x'*Q*x + q'*x;
%funzione gradiente
Df = @(x) 2*Q*x + q;

[K, n] = size(P);

tic
%START
i = 0;
D = Df(x);
y = zeros(n, 1);
fx(1) = f(x);

for k = 1 : K
    Ik = find(P(k,:) == 1); %prendo gli indici diversi da zero (indici appartenenti a Ik)
    Dk = D(Ik); %estraggo le componenti del gradiente
    [~, j] = min(Dk); %calcolo argmin di Dk
    jk = Ik(j); %prendo il j-esimo indice di Ik
    y(jk) = 1; %pongo y(jk) = 1
end
d = y - x; %direzione di decrescita

obj = D'*d; %prodotto scalare tra il gradiente in x e la direzione di decrescita
O(1) = obj;

%LINE SEARCH
if isequal(options1,'LBM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchLBM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end  
elseif isequal(options1,'QBM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchQBM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end 
elseif isequal(options1, 'NM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchNM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end  
else
    alpha = 2/(i + 2);
end

%aggiorno il vettore x
x = x + alpha * d;

%valuto la funzione nel nuovo punto
fx(2) = f(x);

%%%%%%%%%%%%%%%%%MOMENTUM%%%%%%%%%%%%%%%%%%%%
if isequal(options2, 'yes')
    alpha_old = alpha;
    d_old = y - x;
    %d_old = d;
end
%%%%%%%%%%%%%%%%%MOMENTUM%%%%%%%%%%%%%%%%%%%%

i = i + 1;

while (obj < - eps)
    D = Df(x);
    y = zeros(n, 1);
    for k = 1 : K
        Ik = find(P(k,:) == 1);
        Dk = D(Ik);
        [~, j] = min(Dk);
        jk = Ik(j);
        y(jk) = 1;
    end
    d = y - x;
    obj = D'*d;
    O(i) = obj;
    if isequal(options1,'LBM')
        alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
        if (alphaStart <= 1)
            alpha = LineSearchLBM(Q, q, x, d, alphaStart, eps_ls);
        else
            alpha = 1;
        end  
    elseif isequal(options1,'QBM')
        alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
        if (alphaStart <= 1)
            alpha = LineSearchQBM(Q, q, x, d, alphaStart, eps_ls);
        else
            alpha = 1;
        end 
    elseif isequal(options1, 'NM')
        alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
        if (alphaStart <= 1)
            alpha = LineSearchNM(Q, q, x, d, alphaStart, eps_ls);
        else
            alpha = 1;
        end  
    else
        alpha = 2/(i + 2);
    end
    if isequal(options2, 'no')
        x = x + alpha * d;
        fx(i+1) = f(x);
    end
    %%%%%%%%%%%%%%%%%MOMENTUM%%%%%%%%%%%%%%%%%%%%
    if isequal(options2, 'yes')
        %Si impone che il nuovo punto sia all'interno del triangolo avente
        %vertici x, d_old, d_new
        par_momentum = min([alpha_old, 1 - alpha_old, alpha, 1 - alpha]);
        momentum = par_momentum * d_old;
        x = x + alpha * d + momentum;
        fx(i+1) = f(x);
        alpha_old = alpha;
        d_old = y - x;
        %d_old = d;
    end
    %%%%%%%%%%%%%%%%%MOMENTUM%%%%%%%%%%%%%%%%%%%%
    i = i + 1;
end

elapsedTime = toc;
x_min = x;
f_min = f(x);
numSteps = i;

figure('Name','Objective Function');
plot(O, 'ro-')
hold on
plot(fx, 'bo-')
hold off


