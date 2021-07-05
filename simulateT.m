% SIMULATION INPUT

load Body.mat P % Read Body data
load Probes.mat buttons % Read sensor locations
Tamb = 20; % Ambient Temperature [degC]
sigma = 5.6704*10^(-8) ; % Stefan-Boltzmann constant [W/(m^2*K^4)]
emis = 0.96; % Emissivity of human skin
shfl = 3.81; % Specific heat capacity floor [MJ/(m^3*K)]
lfl = 14.35; % Thermal conductivity floor [W/(m*K)]
ws = 0; % Windspeed [m/s]
cl = 1.0525; % Characteristic length of the body for the convection (Calculated as the average of the length of the body and the diameter of the body) [cm]
clinv = 1/cl;
mosst = 0.6397; % Ratio of the surface area of the body (calculated using Mossteller's forumla) and the surface area of the cubic model
rekentijd = 48; % Post-mortem period over which the heat exchange is simulated [h]

pp = getInput(P,Tamb,cl,sigma,emis,shfl,lfl,ws);

neighbours = pp.neighbours;
tijdstap = pp.tijdstap;
Para = pp.Para;
Tenv = pp.Tenv;

LossFr = pp.LossFr;
LUTforced = pp.LUTforced;
LUTfree = pp.LUTfree;
airar = pp.airar;
airmats = pp.airmats;
idxF = pp.idxF;
sN = pp.sN;

L = pp.L;
H = pp.H;

% INITIALIZE VIRTUAL PROBES

bodyparts = fieldnames(buttons);
sB = size(bodyparts,1);

Curves = zeros(sB+1,1);
coordinates = zeros(sB,3);
for uu=1:sB
    coordinates(uu,:) = buttons.(bodyparts{uu,1}).coordinates;
    Curves(uu,1) = L(coordinates(uu,1),coordinates(uu,2),coordinates(uu,3),1) ;
end
Curves(sB+1,1) = Para(2) ;

% SIMULATE HEAT EXCHANGE

for Tijd = 1:(rekentijd*3600/tijdstap)
    for gg = sN:-1:1
        i = neighbours(1,gg);
        j = neighbours(2,gg);
        k = neighbours(3,gg);
        airi = airar(gg);
        hcur = 1/H(i,j,k);
        pcur = P(i,j,k);
        lcur = L(i,j,k,1);
        Loss(1) = 1/(hcur + 1/H(i-1,j,k)) * LossFr(pcur) * (lcur - L(i-1,j,k,1));
        Loss(2) = 1/(hcur + 1/H(i+1,j,k)) * LossFr(pcur) * (lcur - L(i+1,j,k,1));
        Loss(3) = 1/(hcur + 1/H(i,j-1,k)) * LossFr(pcur) * (lcur - L(i,j-1,k,1));
        Loss(4) = 1/(hcur + 1/H(i,j+1,k)) * LossFr(pcur) * (lcur - L(i,j+1,k,1));
        Loss(5) = 1/(hcur + 1/H(i,j,k-1)) * LossFr(pcur) * (lcur - L(i,j,k-1,1));
        Loss(6) = 1/(hcur + 1/H(i,j,k+1)) * LossFr(pcur) * (lcur - L(i,j,k+1,1));
        if airi
            airs = airmats(:,gg);
            hrad = emis*sigma*(lcur^2 + Tenv^2)*(lcur+Tenv);  
            lug = (lcur+Tenv)/2-273.15 ;
            if lug<=12.5
                lug = 1 ;
            elseif lug <= 17.5
                lug = 2 ;
            elseif lug <= 22.5
                lug = 3 ;
            elseif lug <=27.5
                lug = 4 ;
            elseif lug <=32.5
                lug = 5 ;
            else
                lug = 6 ;
            end
            Nuforcedlucht = LUTforced(lug,idxF);
            NufreeL = LUTfree(lug);
            Nufreelucht = NufreeL*abs(lcur-Tenv)^0.25 ;
            Nulucht = (Nufreelucht^3+Nuforcedlucht^3)^(1/3) ;
            hconv = clinv*Nulucht*2.6*10^-2 ;
            LossAir = (1/(hcur + 1/(hrad + hconv)) * LossFr(pcur) * (lcur - Tenv))* mosst;
            Loss(airs) = LossAir;
        end
        L(i,j,k,2) = lcur - sum(Loss) ;
    end

    for uu=1:sB
        Curves(uu,Tijd+1) = L(coordinates(uu,1),coordinates(uu,2),coordinates(uu,3),1) ;
    end
    Curves(sB+1,Tijd+1) = Tenv ;
    for k = 1:size(P,3)
        L(:,:,k,1) = L(:,:,k,2) ;
    end

end

% PLOT SIMULATED TEMPERATURES

figure(1)
plot((1:size(Curves,2))./60,Curves(1:sB,:) - 273.15,'LineWidth',1)
legend(bodyparts)
xlabel('PMI (h)')
ylabel('Temperature (\circ C)')