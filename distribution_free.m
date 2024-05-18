clear all
deltaE = 10;
E = 1:deltaE:1000;
R = 8.314;
beta = 20 / 60;
deltat = 0.1;
load('intermls.mat');
target = mlsf(1:80001,1);
u1 = 130; delta1 = 30; u2 = 550; delta2 = 30;
f = (0.45 * exp(-(E-u1).^2./2/delta1/delta1) / delta1 + 0.55 * exp(-(E-u2).^2./2/delta2/delta2) / delta2) / sqrt(2*pi);
n = 2;
for ii = 1:length(E)
    vistar(ii) = rand() ;
    A(ii) = exp(0.051*E(ii) + 7.8);
    dvi(ii) = 0;
    dA(ii) = 0;  
end
vistar = vistar / sum(vistar) * (1/deltaE);
AA = A;
format long
tic
for iteration = 1:100
    for ii = 1:length(E) 
        vi(1:80001,ii) = 0;
        intAi(1:80001,ii) = 0;
        for t = 0.1:deltat:20000
            T = beta * t;
            k = round(t/deltat);
            if k > 1
                dmls(k,ii) = A(ii)*exp(-E(ii)/R/T*1000) * (vistar(ii) - vi(k-1,ii));
                vi(k,ii) = vi(k-1,ii) + dmls(k,ii) * deltat;    
                intAi(k,ii) = intAi(k-1,ii) + A(ii) * exp(-E(ii)/R/T *1000)*deltat;
            else
                dmls(k,ii) = A(ii)*exp(-E(ii)/R/T*1000) * (vistar(ii) - vi(1,ii));
                vi(k,ii) = vi(k,ii) + dmls(k,ii) * deltat;
                intAi(k,ii) = intAi(k,ii) + A(ii) * exp(-E(ii)/R/T *1000)*deltat;
            end           
            if k > 80000
                break;
            end            
        end
    end
    dmlsf = sum(dmls,2) * deltaE;
    for ii = 1:length(E) 
        dvi(ii) = sum(2 * (dmlsf - target).* dmls(:,ii)) * deltat;
        dA(ii) = sum(2 * (dmlsf - target).* dmls(:,ii).*(1-intAi(:,ii))) * deltat;
    end
    vect = [dvi, dA]';
    dev(iteration) = sum((dmlsf - target).^2)*deltat;
    gradient(iteration) = norm(vect(1:length(E))');
    [iteration dev(iteration) gradient(iteration)]
    step = dev(iteration)/(sqrt(vect'*vect) + 1e-6);
    vistar = vistar - step * (vect(1:length(E))')*10;
    lnA = step * (vect(length(E)+1:end)')* 1e5;
    lineA = [ones(size(lnA)); 1:length(lnA)] * [ones(size(lnA)); 1:length(lnA)]';
    lineB = [ones(size(lnA)); 1:length(lnA)]*lnA';
    coef = inv(lineA)*lineB;
    lnnA = coef(1) + coef(2) * (1:length(lnA));   
   if mod(iteration,50) == 1 && iteration < 2000
        vistar = sgolayfilt(vistar, 3, 5); 
    end
    if iteration > 200
        A = exp(log(A) -lnnA * 10); 
    end
    vistar(vistar < 0) = 0;
    vistar = vistar / sum(vistar) * (1/deltaE);
end
toc