% ----------------------------------------------------------------------- %
%   Author:  Stefano Aretusini                                            %
%   Date:    06/06/2020                                                   %
%   E-mail:  ste.aretu (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %

function Output=last_implicitmodel_diffusivity(Input)

% Input arguments & logics
a=[isfield(Input,'Time'),isfield(Input,'shear1'),isfield(Input,'vel')];

if any(a==0)
    error('Please pr ovide input Time, shear or velocity into the input struct')
else
    Time=Input.Time;
    shear1=Input.shear1;
    vel=Input.vel;
end

if isfield(Input,'phi')
    phi=Input.phi;
else
    phi=ones(size(Time,1),size(Time,2))*0.05;
    disp(['default constant porosity: ',num2str(phi(end)),' %'])
end

if isfield(Input,'totalvolume')
totalvolume=Input.totalvolume; %this is basically not used anymore, check!
else
    totalvolume=ones(size(Time,1),size(Time,2))*(pi*((50/2000)^2)-pi*((30/2000)^2))*0.002;
    disp(['default constant sample volume: ',num2str(totalvolume(end)),' m^3'])
end

% Size of t domain
if isfield(Input,'tsize')
    tsize=Input.tsize;
else
    tsize=Time(end);
    disp(['default size of time domain: ',num2str(tsize),' s'])
end

if isfield(Input,'tnum')
    tnum=Input.tnum;
else
    tnum=101;
    disp(['default number of nodes in time: ',num2str(tnum)])
end

if isfield(Input,'exptime')
exptime=Input.exptime; % this is basically not used anymore, check!
else
    exptime=Time(end);
    disp(['default total time of the experiment: ',num2str(exptime),' s'])
end
% Size of x domain
if isfield(Input,'h0')
    h0=Input.h0;
else
    h0=0.01;
    disp(['default thickness of gouge domain: ',num2str(h0),' m'])
end

if isfield(Input,'xnum')
    xnum=Input.xnum;
else
    xnum=101;
    disp(['default number of nodes in x: ',num2str(xnum)])
end

% Initial conditions
if isfield(Input,'initialT')
    initialT=Input.initialT;
else
    initialT=298;
    disp(['default initial temperature: ',num2str(initialT),' K'])
end

if isfield(Input,'initialP')
    initialP=Input.initialP;
else
    initialP=100000;
    disp(['default initial pressure: ',num2str(initialP),' Pa'])
end

% Boundary conditions
if and(isfield(Input,'Pus_boundary_option'),isfield(Input,'p0us'))
    Pus_boundary_option=Input.Pus_boundary_option;
    p0us=Input.p0us;
else
    Pus_boundary_option='constant value';
    p0us=100000;
    disp(['default ',Pus_boundary_option,' upstream pressure: ',num2str(p0us),' Pa'])
end

if and(isfield(Input,'Pds_boundary_option'),isfield(Input,'p0ds'))
    Pds_boundary_option=Input.Pds_boundary_option;
    p0ds=Input.p0ds;
else
    Pds_boundary_option='constant value';
    p0ds=100000;
    disp(['default ',Pds_boundary_option,' downstream pressure: ',num2str(p0ds),' Pa'])
end

if and(isfield(Input,'Tus_boundary_option'),isfield(Input,'t0us'))
    Tus_boundary_option=Input.Tus_boundary_option;
    t0us=Input.t0us;
else
    Tus_boundary_option='constant value';
    t0us=298;
    disp(['default ',Tus_boundary_option,' upstream temperature: ',num2str(t0us),' K'])
end

if and(isfield(Input,'Tds_boundary_option'),isfield(Input,'t0ds'))
    Tds_boundary_option=Input.Tds_boundary_option;
    t0ds=Input.t0ds;
else
    Tds_boundary_option='constant value';
    t0ds=298;
    disp(['default ',Tds_boundary_option,' downstream temperature: ',num2str(t0ds),' K'])
end

% Thermal and hydraulic properties (solids)
a=[isfield(Input,'K'),isfield(Input,'lambdac'),isfield(Input,'thc'),isfield(Input,'rho'),isfield(Input,'cp')];

if any(a==0)
    K=1e-15;
    lambdac=2e-4;
    thc=1.5;
    rho=2600;
    cp=800;
    
    disp(['default permeability: ',num2str(K),' m^2'])
    disp(['default thermal expansivity: ',num2str(lambdac),' 1/K'])
    disp(['default thermal conductivity: ',num2str(thc),' W/(m*K)'])
    disp(['default density: ',num2str(rho),' kg/m^3'])
    disp(['default heat capacity: ',num2str(cp),' J/(kg*K)'])
    
else
    K=Input.K;
    lambdac=Input.lambdac;
    thc=Input.thc;
    rho=Input.rho;
    cp=Input.cp;
end

% Thermal and hydraulic properties (fluid)

if isfield(Input,'fluidname')
    fluidname=Input.fluidname;
    
    %used to query the tables for properties this needs "table"
    % this should go back as it was (ask for separate input values)
    % optional: include XSteam
else
    fluidname='water';
    disp(['default fluid: ',fluidname,' (room PT conditions)'])
end

% Width of localized area (has to be lower than h0)
if isfield(Input,'width')
    width=Input.width;
else
    width=1e-4;
    disp(['default thickness of localized heat source: ',num2str(width),' m'])
end

% Thermalpressurization active or not
if isfield(Input,'thermalpress_option')
    thermalpress_option=Input.thermalpress_option;
else
    thermalpress_option=0;
    disp('default: thermal pressurization not active')
    
end

% Thermochemical active if reactionname~='none'
if isfield(Input,'reactionname')
    reactionname=Input.reactionname;
else
    reactionname='none';
    disp('default: thermal decomposition not active')
end

% Shear compaction contribution calculated from:
% Porosity
if isfield(Input,'porosity_option')
    porosity_option=Input.porosity_option;
else
    porosity_option=0;
    disp('default: shear compaction from porosity change not active')
end

% Pump volume change
if and(isfield(Input,'pump_option'),isfield(Input,'P1'))
    pump_option=Input.pump_option;
    P1=Input.P1;
else
    pump_option=0;
    P1=0;
    disp('default: shear compaction from volumomemetry not active')
end

% a trick to limit max pf equal to sigma_n
if and(isfield(Input,'damp_option'),isfield(Input,'sntot'))
    damp_option=Input.damp_option;
    sntot=Input.sntot;
else
    damp_option=0;
    sntot=0;
    disp('default: no tricks ;)')
end

% Option to add additional nodes around gouge layer with mottcorp porous
% plates properties
if and(isfield(Input,'ForcingBlock'),isfield(Input,'FBmaterialname'))
    ForcingBlock=Input.ForcingBlock;
    FBmaterialname=Input.FBmaterialname;
    
else
    ForcingBlock=0;
    FBmaterialname=0;
    disp('default: no forcing block around sample')
end

% T-dependent rho and cp
if and(isfield(Input, 'variableTthermalprop'),isfield(Input,'materialname'))
    variableTthermalprop=Input.variableTthermalprop;
    materialname=Input.materialname; % required by functions "heatcapacity" and "thermalconductivity" called if variableTthermalprop==1
else
    variableTthermalprop=0;
    materialname='';
    disp('default: no temperature dependent thermal properties of solids')
end
% NB: conservation of T and P flux is included if metalFB==1 or variableTthermalprop==1

if isfield(Input,'tables')
    tables=Input.tables;
else
    
    sheet={'fluidproperties','properties','heatcapacity','thermalconductivity','reactions','flowlaws'};
    fname = {'fluidproperties','properties','heatcapacity','thermalconductivity','reactions','flowlaws'};
    for k=1:size(fname,2)
        opts = detectImportOptions('rock_properties.xlsx','Sheet',sheet{k});
        
        for i=1:size(opts.VariableTypes,2)
            if opts.VariableTypes{1,i}=="char"
                opts.VariableTypes{1,i}='categorical';
            else
            end
        end
        
        opts.Sheet=sheet{k};
        opts.DataRange='A3';
        opts.VariableNamesRange='1:1';
        opts.VariableUnitsRange='2:2';
        
        
        tables.(fname{k}) = readtable('rock_properties.xlsx',opts);
        
    end
    
    disp('loaded table from rock_properties.xlsx')
end

if isfield(Input,'equivalent_radius')
    equivalent_radius=Input.equivalent_radius;
else
    ri=30/2000;
    ro=50/2000;
    equivalent_radius=2/3*(ri^2+ri*ro+ro^2)/(ri+ro);
    disp('default equivalent radius for 30/50 sample holder')
end

if isfield(Input,'radius')
    radius=Input.radius;
else
    ri=30/2000;
    ro=50/2000;
    radius=2/3*(ri^2+ri*ro+ro^2)/(ri+ro);
    disp('default radius equal to equivalent radius for 30/50 sample holder')
end

if isfield(Input,'StrainRateFunction')
    StrainRateFunction=Input.StrainRateFunction;
else
    StrainRateFunction='dirac';
    disp('default strain rate function is a regularized dirac delta')
end

if isfield(Input,'ArAn')
    ArAn=Input.ArAn;
else
    ArAn=1;
    disp('default Ar/An is 1')
end

if isfield(Input,'sym')
    sym=Input.sym;
else
    sym=0.5;
    disp('default sym is 0.5')
end

% NB: the single functions open the tables if not included as input here only for fluids this does not happen.

%% NODES PROPERTIES AND TIMESTEP

xsize=h0;
xstp=xsize/(xnum-1);

% re-parametrize the x-domain
% |--downstreamFB--|-----xsize-----|--upstreamFB--|
% in xsize: |-----!-----| heat source at center: sym=1/2
% in xsize: |---------!-| heat source at center: sym=1/10

% xsize domain
L1=-(xsize-sym*xsize);
L2=sym*xsize;

x=L1:xstp:L2;

if ForcingBlock==1
    
    %set this reasonably large to avoid that the boundary conditions affect
    %the modelled temperature (it was 0.02)
    xsizeus=0.1;
    xsizeds=0.1;
    
    xFBu=L1-xsizeus:xstp:L1;
    xFBd=L2:xstp:L2+xsizeds;
    
    %append the forcing block and update x, xnum
    x=[xFBu x xFBd];
    xnum=size(x,2);
end

dt=tsize/(tnum-1);
texp=exptime/dt;


%% INITIAL CONDITIONS AND SOURCE

if initialT==0, t0=298; else, t0=initialT; end
if initialP==0, p0=1e5; else, p0=initialP; end

p0imp=ones(xnum,1)*p0;
t0imp=ones(xnum,1)*t0;

if pump_option==1
    disp(P1)
else
    P1=0;
end

%% INITIALIZATION OF VECTORS

porosity_exp=zeros(tnum,1);

damp=zeros(xnum,tnum);
damp2=zeros(xnum,tnum);

QP_shcomp=zeros(tnum,1);
QP_tp=zeros(xnum,1);
% QT1=zeros(tnum,1);

app.model.Time=zeros(tnum,1);
app.model.phi=zeros(tnum,1);
app.model.betac=zeros(tnum,1);
app.model.kappa_hy=zeros(tnum,1);
app.model.Lambda=zeros(tnum,1);

app.model.pressure=zeros(xnum,tnum);
app.model.temperature=zeros(xnum,tnum);

app.model.QP=zeros(xnum,tnum);
app.model.QP_nochem=zeros(xnum,tnum);
app.model.QP_shcomp=zeros(xnum,tnum);
app.model.QP_tp=zeros(xnum,tnum);
app.model.damp=zeros(xnum,tnum);
% app.model.Psource=zeros(xnum,tnum);

app.model.QT=zeros(xnum,tnum);
app.model.QT_fh=zeros(xnum,tnum);
% app.model.QT1=zeros(tnum,1);
% app.model.Tsource=zeros(xnum,tnum);

%% COUNTERS
timesum=0; % Elapsed time

%% INTERPOLATE INPUT VECTORS TO MODEL TIME VECTOR
app.model.Time=0:dt:tsize;

shear1mod = interp1(Time,shear1,app.model.Time,'linear');
velmod = interp1(Time,vel,app.model.Time,'linear');
phimod = interp1(Time,phi,app.model.Time,'linear');

shear1mod(isnan(shear1mod))=0; shear1mod=shear1mod/ArAn;
velmod(isnan(velmod))=0;
phimod(isnan(phimod))=0;

%angular velocity
ang_vel_mod=velmod/(2*pi*equivalent_radius);

%convert to new radius
new_velmod=2*pi*radius*ang_vel_mod;
velmod=new_velmod;

app.model.shear1=shear1mod;
app.model.vel=velmod;
app.model.phi=phimod;
app.model.ang_vel=ang_vel_mod;


%% TIME LOOP
for t=1:1:tnum
    
    %% PHISYCAL PROPERTIES
    
    % THERMAL PROPERTIES (optional: T-dependent)
    
    % rho, cp, thc are defined for the gouge
    % rhoFB, cpFB, thcFB are defined for the Forcing Block
    
    if variableTthermalprop==1
        
        cp=heatcapacity(tables,t0imp,materialname);
        thc=thermalconductivity(tables,t0imp,materialname);
        
        rhocp=rho.*cp; % volumetric heat capacity, j/kg/m^3
        kappa=thc./rhocp; % Thermal diffusivity, m^2/s
        
        KT=kappa;
        RCT=rhocp;
        
    else
        rhocp=rho.*cp; % volumetric heat capacity, j/kg/m^3
        kappa=thc./rhocp; % Thermal diffusivity, m^2/s
        
    end
    
    % thermal properties of the forcing block
    if ForcingBlock==1
        
        switch FBmaterialname
            case 'mottcorp frit'
                % data from mottcorp grade 0.5 (grade 5 was preloaded)
                rhoFB=6162;
                cpFB=475;
                thcFB=8.6;
                
            case 'carrara rock'
                rhoFB=2860;
                cpFB=830;
                thcFB=1.8;
                
            case 'WC' %http://www.allaboutcementedcarbide.com/
                rhoFB=15700; %pure endmember
                cpFB=350; %150-350
                thcFB=80; %80 fine, 92 medium, 105 coarse
                
            case 'steel' %aisi304, http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MQ304A
                rhoFB=8000;
                cpFB=500;
                thcFB=16.2;
                
            case 'Ti6Al4V' %https://www.azom.com/properties.aspx?ArticleID=1547
                rhoFB=4512;
                cpFB=570;
                thcFB=7.3;
        end
        
        rhocpFB=rhoFB*cpFB;
        kappaFB=thcFB/(rhoFB*cpFB);
        
        % conductivity changing with x
        KT=ones(xnum,1)*kappaFB;
        RCT=ones(xnum,1)*rhocpFB;
        
        %set conductivity and heat capacity in gouge domain
        for i=1:xnum
            if variableTthermalprop==1
                if and(x(i)>L1,x(i)<L2)
                    KT(i)=kappa(i);
                    RCT(i)=rhocp(i);
                end
            else
                if and(x(i)>L1,x(i)<L2)
                    KT(i)=kappa;
                    RCT(i)=rhocp;
                end
            end
        end
    end
    
    %HYDRAULIC PROPERTIES (time dependent)
    row=tables.fluidproperties.fluidname==fluidname;
    betaf=tables.fluidproperties(row,'betaf'); app.model.betaf=betaf{1,1};
    visc=tables.fluidproperties(row,'visc'); app.model.visc=visc{1,1};
    lambdaf=tables.fluidproperties(row,'lambdaf'); app.model.lambdaf=lambdaf{1,1};
    
    % qui forse si può usare direttamente phimod e non phi!
    if t<texp
        if timesum==0
            porosity=phi(1);
        else
            porosity=phimod(t);
        end
    else
        porosity=phi(end);
    end
    
    beta_c=porosity*app.model.betaf; %1/Pa
    kappa_hy=K/app.model.visc/beta_c; % Hydraulic diffusivity, m^2/s
    Lambda=(app.model.lambdaf-(lambdac))/(beta_c+app.model.betaf); % Pa/K, rice et al., 2006, table b
    
    % hydraulic properties of the forcing block
    
    if ForcingBlock==1
        
        switch FBmaterialname
            case 'mottcorp frit'
                % mottcorp grade 0.5 (grade 5 was preloaded)
                porosityFB=0.21;
                permeabilityFB=1.25e-12; %for gas: permeabilityFB=9.8e-13;
                lambdacFB=0;
                
            case 'carrara rock'
                porosityFB=0.05;
                permeabilityFB=1e-19;
                lambdacFB=0.0002;
                
            case 'WC'
                porosityFB=0.01;
                permeabilityFB=1e-21;
                lambdacFB=4.4e-6; %https://nvlpubs.nist.gov/nistpubs/jres/18/jresv18n1p47_A1b.pdf
                
            case 'steel'
                porosityFB=0.01;
                permeabilityFB=1e-21;
                lambdacFB=17.3e-6; %http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MQ304A
                
            case 'Ti6Al4V'
                porosityFB=0.01;
                permeabilityFB=1e-21;
                lambdacFB=9.1e-6; %https://www.azom.com/properties.aspx?ArticleID=1547
        end
        
        beta_cFB=porosityFB*app.model.betaf;
        kappa_hyFB=permeabilityFB/app.model.visc/beta_cFB;
        LambdaFB=(app.model.lambdaf-lambdacFB)/(beta_cFB+app.model.betaf);
        
        %         conductivity changing with x
        KP=ones(xnum,1)*kappa_hyFB;
        LAMBDAP=ones(xnum,1)*LambdaFB;
        
        %         set conductivity and heat capacity in gouge domain
        for i=1:xnum
            if and(x(i)>L1,x(i)<L2)
                KP(i)=kappa_hy;
                LAMBDAP(i)=Lambda;
            end
        end
    end
    
    %     KP=ones(xnum,1)*kappa_hy;
    
    %% TEMPERATURE SOURCE
    
    %% temperature source by shear heating - gouge layer
    
    switch StrainRateFunction
        case 'gauss'
            Fnew=gauss(x,width)*velmod(t);
        case 'dirac'
            Fnew=delta_reg(x)*velmod(t);
    end
    
    %BULK HEAT SOURCES ONLY! NO FLASH HEATING STUFF!
    %heat source versus t
    
    if or(ForcingBlock==1,variableTthermalprop==1)
        
        %source in x,t (frictional heating)
        QT_fh=abs(1e6*shear1mod(t))*(1./RCT.*Fnew');
        
        if reactionname~="none"
                        
            QT_chem=TemperatureSource(tables,dt,t0imp,reactionname,porosity,1)./RCT; 
            alpha=reacted(tables,timesum,t0imp,reactionname);
            dalphadt=reactedrate(tables,timesum,t0imp,reactionname);
            
            %set QT_chem = 0 in FB nodes.
            [~,idxu]=intersect(x,xFBu,'stable');
            [~,idxd]=intersect(x,xFBd,'stable');
            QT_chem(idxu)=0;
            QT_chem(idxd)=0;
            QT_chem(x==L1)=0;
            QT_chem(x==L2)=0;

            QT=QT_fh+QT_chem;
            
        else
            QT_chem=zeros(size(QT_fh,1),size(QT_fh,2));
            alpha=zeros(size(QT_fh,1),size(QT_fh,2));
            dalphadt=zeros(size(QT_fh,1),size(QT_fh,2));
            
            QT=QT_fh;
        end
        
        
    else       
        
        %source in x,t (frictional heating)
        QT_fh=abs(1e6*shear1mod(t))*1/rhocp*Fnew';
        
        if reactionname~="none"
                        
            QT_chem=TemperatureSource(tables,dt,t0imp,reactionname,porosity,1)/rhocp;
            alpha=reacted(tables,timesum,t0imp,reactionname);
            dalphadt=reactedrate(tables,timesum,t0imp,reactionname);
            
            QT=QT_fh+QT_chem;
            
        else
            QT_chem=zeros(size(QT_fh,1),size(QT_fh,2));
            alpha=zeros(size(QT_fh,1),size(QT_fh,2));
            dalphadt=zeros(size(QT_fh,1),size(QT_fh,2));
            
            QT=QT_fh;
        end
        
    end
    
    %% TEMPERATURE SOLVER WITH OPTIONS
    
    
    if or(ForcingBlock==1,variableTthermalprop==1)
        % left side
        i=2:1:xnum-1;
        
        LT=sparse(i,i-1,-(KT(i-1,1)+KT(i,1))/2/xstp^2,xnum,xnum)+...
            sparse (i,i,1/dt+(KT(i-1,1)+KT(i,1))/2/xstp^2+(KT(i+1,1)+KT(i,1))/2/xstp^2,xnum,xnum)+...
            sparse(i,i+1,-(KT(i+1,1)+KT(i,1))/2/xstp^2,xnum,xnum);
        % right side
        RT=zeros(xnum,1);
        for i=2:1:xnum-1
            RT(i,1)=t0imp(i)/dt+QT(i);
        end
    else %constant properties across x
        %nodes in the middle
        i=2:1:xnum-1;
        
        %left side
        LT=sparse(i,i-1,-kappa/xstp^2,xnum,xnum)+ ...
            sparse(i,i,1/dt+2*kappa/xstp^2,xnum,xnum)+ ...
            sparse(i,i+1,-kappa/xstp^2,xnum,xnum);
        
        %right side
        RT=zeros(xnum,1);
        for i=2:1:xnum-1
            RT(i,1)=t0imp(i)/dt+QT(i);
        end
        
    end
    
    %first node upstream
    switch Tus_boundary_option
        case 'constant value'
            LT=LT+sparse(1,1,1,xnum,xnum);
            RT(1,1)=t0us;
        case 'no flow'
            LT=LT+sparse(1,1,1,xnum,xnum)+sparse(1,1+1,-1,xnum,xnum);
            RT(1,1)=0;
    end
    
    %last node downstream
    switch Tds_boundary_option
        case 'constant value'
            LT=LT+sparse(xnum,xnum,1,xnum,xnum);
            RT(xnum,1)=t0ds;
        case 'no flow'
            LT=LT+sparse(xnum,xnum,1,xnum,xnum)+sparse(xnum,xnum-1,-1,xnum,xnum);
            RT(xnum,1)=0;
    end
    
    %solver
    t1imp=LT\RT;
    
    
    %% PRESSURE SOURCES
    
    %dampening the compaction depending on overpressures
    if damp_option==1
        
        damp(:,t)=(p0imp-p0)./(sntot-p0imp);
        damp2(:,t)=p0imp/sntot;
        
        
        for q=1:size(damp,1)
            
            if abs(damp(q,t))>= 1
                damp(q,t)=1;
            else
            end
            if abs(damp2(q,t))>=1
                damp2(q,t)=1;
            else
            end
            
        end
    else
        damp=zeros(xnum,t);
        damp2=zeros(xnum,t);
    end
    
    %scalar source terms (with varying timestep)
    
    if timesum==0
        P1_phi=0;
    else
        P1_phi=-(porosity-porosity_exp(t-1))/dt;
    end
    
    P1_pump=1e-6*P1/totalvolume(end)/texp;
    
    if t<texp
        
        %shear compaction source term, Pa/s
        
        
        if and(pump_option==1, porosity_option==0)
            QP_shcomp(t)=P1_pump/beta_c;
        elseif and(pump_option==0,porosity_option==1)
            QP_shcomp(t)=P1_phi/beta_c;
        else
            QP_shcomp(t)=0;
        end
        
        %thermal pressurization source term, K/s*Pa/K = Pa/s
        if thermalpress_option==1
            if or(ForcingBlock==1,variableTthermalprop==1)
                QP_tp(:,t)=(t1imp-t0imp)/dt.*LAMBDAP;
            else
                QP_tp(:,t)=(t1imp-t0imp)/dt*Lambda;
            end
        else
            QP_tp(:,t)=zeros(size(Fnew,1),1);
        end
        
    else
        QP_shcomp(t)=0;
        
        QP_tp(:,t)=zeros(size(Fnew,1),1);
        
    end
    
    QP_nochem=(Fnew'*QP_shcomp(t)+QP_tp(:,t)).*(1-damp2(:,t));
    
    
    if reactionname~="none"

        QP_chem=PressureSource(tables,dt,t0imp,reactionname,porosity,1)/beta_c;
        QP=QP_nochem+QP_chem;
    else
        QP_chem=zeros(size(QP_nochem,1),size(QP_nochem,2));
        QP=QP_nochem;
    end
    
    %% PRESSURE SOLVER WITH OPTIONS
    
    if ForcingBlock==1
        % left side
        i=2:1:xnum-1;
        
        LP=sparse(i,i-1,-(KP(i-1,1)+KP(i,1))/2/xstp^2,xnum,xnum)+...
            sparse (i,i,1/dt+(KP(i-1,1)+KP(i,1))/2/xstp^2+(KP(i+1,1)+KP(i,1))/2/xstp^2,xnum,xnum)+...
            sparse(i,i+1,-(KP(i+1,1)+KP(i,1))/2/xstp^2,xnum,xnum);
        
        % right side
        RP=zeros(xnum,1);
        for i=2:1:xnum-1
            RP(i,1)=p0imp(i)/dt+QP(i);
        end
        
    else %constant properties across x
        i=2:1:xnum-1;
        LP=sparse(i,i-1,-kappa_hy/xstp^2,xnum,xnum)+...
            sparse(i,i,1/dt+2*kappa_hy/xstp^2,xnum,xnum)+...
            sparse(i,i+1,-kappa_hy/xstp^2,xnum,xnum); %diffusion
        
        RP=zeros(xnum,1);
        for i=2:1:xnum-1
            RP(i,1)=p0imp(i)/dt+QP(i);
        end
        
    end
    
    switch Pus_boundary_option
        case 'constant value'
            LP=LP+sparse(1,1,1,xnum,xnum);
            RP(1,1)=p0us;
        case 'no flow'
            LP=LP+sparse(1,1,1,xnum,xnum)+sparse(1,1+1,-1,xnum,xnum);
            RP(1,1)=0;
    end
    
    switch Pds_boundary_option
        case 'constant value'
            LP=LP+sparse(xnum,xnum,1,xnum,xnum);
            RP(xnum,1)=p0ds;
        case 'no flow'
            LP=LP+sparse(xnum,xnum,1,xnum,xnum)+sparse(xnum,xnum-1,-1,xnum,xnum);
            RP(xnum,1)=0;
    end
    
    p1imp=LP\RP;
    
    %% SAVING DATA FOR POSTPROCESSING
    %% MODEL OUTPUTS DELIVERED TO app.model
    porosity_exp(t)=porosity;
    app.model.phi(t)=porosity;
    
    app.model.betac(t,1)=beta_c;
    app.model.kappa_hy(t,1)=kappa_hy;
    app.model.Lambda(t,1)=Lambda;
    
    if or(ForcingBlock==1,variableTthermalprop==1)
        app.model.RCT(:,t)=RCT;
        app.model.KT(:,t)=KT;
    end
    
    app.model.pressure(:,t)=p1imp;
    app.model.temperature(:,t)=t1imp;
    
    app.model.QP(:,t)=QP;
    app.model.QP_nochem(:,t)=QP_nochem;
    app.model.QP_shcomp(:,t)=QP_shcomp(t)*Fnew';
    app.model.QP_chem(:,t)=QP_chem;
    
    app.model.QT(:,t)=QT;
    app.model.QT_fh(:,t)=QT_fh;
    app.model.Fnew(:,t)=Fnew;
    app.model.QT_chem(:,t)=QT_chem;
    
    app.model.alpha(:,t)=alpha;
    app.model.dalphadt(:,t)=dalphadt;
    
    %% COUNTERS; ET CETERA
    timesum=timesum+dt;
    app.statusGauge.Value=round(timesum/tsize)*100;
    
    % Reassign pressure profiles for the next step
    p0imp=p1imp;
    t0imp=t1imp;
    
end

%% MODEL OUTPUTS DELIVERED TO app.model
%             app.model.QP_shcomp=QP_shcomp;
app.model.QP_tp=QP_tp;
app.model.damp=damp;
% app.model.QT1=QT1;

% single values
app.model.x=x;
app.model.t=t;
app.model.xsize=xsize;
app.model.tsize=tsize;
app.model.xnum=xnum;
app.model.tnum=tnum;
app.model.p0=p0;
app.model.t0=t0;
app.model.dt=dt;
app.model.dx=xstp;
app.model.sym=sym;

%%GET SUBSETS OF QP_shcomp, QP, P, QT, T solutions

% at center:
[~,i]=min(abs(x-0));

app.model.QP_shcompmid=app.model.QP_shcomp(i,:);
app.model.QP_tpmid=app.model.QP_tp(i,:);
app.model.QP_mid=app.model.QP(i,:);
app.model.P_mid=app.model.pressure(i,:);

app.model.QT_mid=app.model.QT(i,:);
app.model.QT_fh_mid=app.model.QT_fh(i,:);
app.model.QT_chem_mid=app.model.QT_chem(i,:);
app.model.T_mid=app.model.temperature(i,:);

if ForcingBlock==1
    %PRESSURE: position at L2 (upstream) and L1 (downstream)
    [~,i]=min(abs(x-(L2)));
    app.model.Pus=app.model.pressure(i,:);
    
    [~,i]=min(abs(x-(L1)));
    app.model.Pds=app.model.pressure(i,:);
    
    %TEMPERATURE: position at L2 (upstream) and L1 (downstream)
    [~,i]=min(abs(x-(L2)));
    app.model.Tus=app.model.temperature(i,:);
    
    [~,i]=min(abs(x-(L1)));
    app.model.Tds=app.model.temperature(i,:);
else
    app.model.Pus=app.model.pressure(2,:);
    app.model.Pds=app.model.pressure(end-1,:);
    
    app.model.Tus=app.model.temperature(2,:);
    app.model.Tds=app.model.temperature(end-1,:);
end

Output=app.model;

end