% Function: Multi-cell heat microstructure optimization with compatibility
% Author: Tianjie Chen (ctj19990927@mail.ustc.edu.cn)
% Version: 2024-10-13
% V3.0: Adapting new opt model to ensure functionally graded heat metas
% Now can manually assign left and right heat tensors
% Give up using robust TO since there is no confirmed relation between robust and heat targets

% Example  [xPhys_cat, c, Q] = topX_Dual_heat_V2(3,1,100,0.16,31,0.5,2,1,0,1,2000,0.0001,1,1,[10,10],[20,20]);

function [xPhys_cat, c, Q, c_hist,cons_hist] = topX_Dual_heat_V3(CellCount,CELL_TYPE,nel,Kmin,Kmax,w1,GLOBAL,LOCAL,FilterDom,nloop,epsilon,IC,assign,left,right)


%% INPUTS
% CELL-TYPE 1: ISOTROPIC; 2: ORTHOTROPIC; OTHERS: ANISOTROPIC
% OBJECTIVE 1: Bulk Modulus; 2: Shear Modulus; 3: Negative Poisson
% CONSTRAINT TOGGLES (GLOBAL, LOCAL, ISO): enabled(1), disabled(0)
% GLOBAL SCALE FILTERING TOGGLE (FilterDom): enabled(1), disabled(0)

%% NO ROBUST TO SINCE HEAT MATERIALS HAVE NO TARGETS LIKE BULK MODULUS
eta = 0.5;% Here it's no longer parameter for robust TO, but just for projection control
beta = 1;
BetaUpdFreq = 200;
BetaUpdFactor = 2;

%% LOCAL VOLUME CONSTRAINRT PARAMETERS
p = 16;                         % pNorm
%r_hat = 10;                     % pNorm radius: sets maximum length scale (solid material)
r_hat = 6;                     % pNorm radius: sets maximum length scale (solid material)
v = 0.95;

%% MATERIAL PARAMETERS

%% MISC PARAMETERS
penal = 5;
rmin = 3;
%rmin = 2.2;
sym = 4;                                % 0: None 1: 4-fold diagonal 2: 4-fold normal 3: Chiral 4: Cubic
Influence = CellCount;                  % Order of influence; 1: Individual, 2: Adjacent cell compatibility
%Influence = 2;                  % Order of influence; 1: Individual, 2: Adjacent cell compatibility
%IC = 1;
%IC = 2;

if w1 < 0
    error('Illegal weight input!');
end

if w1 == 0
    Influence = 1;
end

nelx = nel;
nely = nel;

CaseCount = Influence*CellCount - Influence*(Influence-1)/2;
kalpa=cell(CellCount,1);
if CELL_TYPE==1
    sym = 4;
    if assign
        kalpafrac_a=linspace(left(1),right(1),CellCount);
        kalpafrac_b=linspace(left(2),right(2),CellCount);
        for i=1:CellCount
            kalpa{i}=[kalpafrac_a(i);kalpafrac_b(i)];
        end
    else
        kalpafrac = linspace(Kmin,Kmax,CellCount+2);
        for i=1:CellCount
            kalpa{i}=[kalpafrac(i+1);kalpafrac(i+1)];
        end
    end
elseif CELL_TYPE==2
    sym = 2;
    if assign
        kalpafrac_a=linspace(left(1),right(1),CellCount);
        kalpafrac_b=linspace(left(2),right(2),CellCount);
        for i=1:CellCount
            kalpa{i}=[kalpafrac_a(i);kalpafrac_b(i)];
        end
    else
        kalpafrac = linspace(Kmin,Kmax,CellCount+2);
        for i=1:CellCount
            kalpa{i}=[kalpafrac(i+1);kalpafrac(i+1)];
        end
    end
else

end

vol_max = linspace(v,v,CellCount);
vol_max_pNorm = (nelx*nely*vol_max.^p).^(1/p);
w = w1;
%w = 100*w;

% used for importing existing design
x_imp = 0;
imp = 0;
% x_imp = imread("experiment5_multiscale_minthermalcompliance/compare_Xcells(10,[5,25]).png");
% x_imp=1-x_imp;
% imp = 1;

%% PREPARE FINITE ELEMENT ANALYSIS
KE = [ 2/3 -1/6 -1/3 -1/6
       -1/6 2/3 -1/6 -1/3
       -1/3 -1/6 2/3 -1/6
       -1/6 -1/3 -1/6 2/3];%heat conductivity matrix of basic cell

nodenrs = cell(2,1);
edofMat = cell(2,1);
iK = cell(2,1);
jK = cell(2,1);
[nodenrs{1}, edofMat{1}, iK{1}, jK{1}] = PrepareFEA(nelx,nely);
[nodenrs{2}, edofMat{2}, iK{2}, jK{2}] = PrepareFEA(2*nelx,nely);


%% PREPARE PDE FILTER FOR LOCAL VOLUME CONSTRAINT
edofVecF = reshape(nodenrs{1}(1:end-1,1:end-1),nelx*nely,1);
edofMatF = repmat(edofVecF,1,4)+repmat([0 nely+[1:2] 1],nelx*nely,1);
iKF = reshape(kron(edofMatF,ones(4,1))',16*nelx*nely,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nelx*nely,1);
iTF = reshape(edofMatF,4*nelx*nely,1);
jTF = reshape(repmat([1:nelx*nely],4,1)',4*nelx*nely,1);
sTF = repmat(1/4,4*nelx*nely,1);
TF = sparse(iTF,jTF,sTF);

Rmin = r_hat/2/sqrt(3);
KEF = Rmin^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
             [4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36;
sKF = reshape(KEF(:)*ones(1,nelx*nely),16*nelx*nely,1);
KF = sparse(iKF,jKF,sKF);
LF = chol(KF,'lower');

%% PERIODIC BOUNDARY CONDITIONS
T_i = cell(2,1);
Tfixed = cell(2,1);
d1 = cell(2,1);
d2 = cell(2,1);
d3 = cell(2,1);
d4 = cell(2,1);
wfixed = cell(2,1);
[T_i{1}, Tfixed{1}, d1{1}, d2{1}, d3{1}, d4{1}, wfixed{1}] = ApplyPeriodicBC(nelx,nely,nodenrs{1,1});
[T_i{2}, Tfixed{2}, d1{2}, d2{2}, d3{2}, d4{2}, wfixed{2}] = ApplyPeriodicBC(2*nelx,nely,nodenrs{2,1});

%% INITIALIZE ITERATION
Q_i = cell(CaseCount,1);
dQ_i = cell(CaseCount,1);

x = cell(CaseCount,1);
xTilde = cell(CaseCount,1);
xPhys_i = cell(CaseCount,1);
xold1 = cell(CaseCount,1);
xold2 = cell(CaseCount,1);

for i = 1:CaseCount
    Q_i{i,1} = zeros(2,2);
    dQ_i{i,1} = cell(2,2);
end

volfrac=linspace(0,1,CellCount+2);
volfrac=volfrac(2:CellCount+1);
for i = 1:CellCount
    x_imp_temp = zeros(nely,nelx);
    if imp == 1
        x_imp_temp = x_imp(:,(1+(i-1)*nelx:i*nelx));
    end
    [x{i}, xTilde{i}, xPhys_i{i}, xold1{i}, xold2{i}] = InitializeIteration(nelx,nely,beta,eta,volfrac(i),IC,x_imp_temp,imp);
end

i = CellCount + 1;
l = 1;
while i <= CaseCount
    j = 1;
    k = j + l;
    while k <= CellCount
        xTilde{i} = horzcat(xTilde{j},xTilde{k});
        xPhys_i{i} = horzcat(xPhys_i{j},xPhys_i{k});
        xold1{i} = horzcat(xold1{j},xold1{k});
        xold2{i} = horzcat(xold2{j},xold2{k});

        i  = i + 1;
        j = j + 1;
        k = j + l;
    end
    l = l + 1;
end

loop = 0;
loopbeta = 0;
change = 1;
low = 0;
upp = 0;


%% START ITERATION
% store results
c_hist = zeros(nloop,CaseCount);        % objective
vol_hist = zeros(nloop,CaseCount);      % volume
change_hist = zeros(nloop,1);   % maximum design change
sharp_hist = zeros(nloop,1);    % sharpness
cons_hist = zeros(nloop,2);     % constraints

plotadd=strcat(num2str(loop,'%04d'),'.png');
imwrite(1-horzcat(x{1:CellCount}),strcat('objective_history_demo\mutual,w=1\iteration_fig\',plotadd));
tic;
%% START ITERATION
while (change > 0.0001) && loop<nloop
    loop = loop + 1;
    loopbeta = loopbeta + 1;

    %% PREPARE FILTER
    iH = ones(CellCount*nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for i1 = 1:CellCount*nelx
        for j1 = 1:nely
            e1 = (i1-1)*nely+j1;
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),CellCount*nelx)
                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                    e2 = (i2-1)*nely+j2;
                    k = k+1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                end
            end
        end
    end
    H = sparse(iH,jH,sH);
    Hs = sum(H,2);

    iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1-1)*nely+j1;
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                    e2 = (i2-1)*nely+j2;
                    k = k+1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                end
            end
        end
    end
    H_ind = sparse(iH,jH,sH);
    Hs_ind = sum(H_ind,2);

    %% FE ANALYSIS
    for i = 1:CaseCount
        if i <= CellCount
            T_i{i} = FEAnalysis(nelx,nely,KE,Kmin,xPhys_i{i},penal,Kmax,iK{1},jK{1},d1{1},d2{1},d3{1},d4{1},Tfixed{1},wfixed{1});
            [Q_i{i}, dQ_i{i}] = Homogenization(nelx,nely,T_i{i},edofMat{1},Kmin,Kmax,xPhys_i{i},penal,KE);
        else
            T_i{i} = FEAnalysis(2*nelx,nely,KE,Kmin,xPhys_i{i},penal,Kmax,iK{2},jK{2},d1{2},d2{2},d3{2},d4{2},Tfixed{2},wfixed{2});
            [Q_i{i}, dQ_i{i}] = Homogenization(2*nelx,nely,T_i{i},edofMat{2},Kmin,Kmax,xPhys_i{i},penal,KE);
        end
    end

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    c_i = 0;
    dc_i = zeros(nely,CellCount*nelx);
    constraints = zeros(CellCount,1);
    dconstraints = cell(CellCount,1);
    current_vol = 0;
    current_heatdamps = 0;
    %Objective functions for individual cells: minimize cell volumes
    %Constraint functions for individual cells: target heat tensor tolerance
    for i = 1:CellCount
        c_i_temp = sum(sum(xPhys_i{i}))/(nelx*nely);
        dc_i_temp = ones(nely,nelx)/(nelx*nely);
        %[c_i_temp,dc_i_temp] = Objective(Q_i{i},dQ_i{i},OBJECTIVE,1);
        [constraints(i),dconstraints{i}] = Constraint(Q_i{i},dQ_i{i},kalpa{i},epsilon);
        c_hist(loop,i) = c_i_temp;
        cons_hist(loop,i)=constraints(i);
        
        c_i =  c_i + c_i_temp;
        current_vol = current_vol+c_i_temp;
        dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) + dc_i_temp;
   end
   %current_vol=current_vol/(CellCount);

   %Objectives for intercells: minimize interccell heat damps
    i = CellCount + 1;
    l = 1;
    while i <= CaseCount
        j = 1;
        k = j + l;
        while k <= CellCount

            [c_i_temp,dc_i_temp] = Objective(Q_i{i},dQ_i{i},2,2);

            c_hist(loop,i) = c_i_temp;

            c_i = c_i + w * c_i_temp;
            current_heatdamps = current_heatdamps+c_i_temp;
            tmp_i = + w * dc_i_temp;
            tmp_i_L = tmp_i(1:nely,1:nelx);
            tmp_i_R = tmp_i(1:nely,(nelx+1):end);
            dc_i(1:nely,(((j-1)*nelx)+1):(j*nelx)) = dc_i(1:nely,(((j-1)*nelx)+1):(j*nelx)) + tmp_i_L;
            dc_i(1:nely,(((k-1)*nelx)+1):(k*nelx)) = dc_i(1:nely,(((k-1)*nelx)+1):(k*nelx)) + tmp_i_R;

            i  = i + 1;
            j = j + 1;
            k = j + l;
        end
        l = l + 1;
    end

    %dv = ones(nely,CellCount*nelx);

    x_cat = horzcat(x{1:CellCount});
    xTilde_cat = horzcat(xTilde{1:CellCount});

    x_pde_hat_i = cell(CellCount,1);
    dfdx_pde_i = cell(CellCount,1);
    for i = 1:CellCount
        x_pde_hat_i{i} = (TF'*(LF'\(LF\(TF*xPhys_i{i}(:)))));
        dfdx_pde_i{i} = (sum(x_pde_hat_i{i}.^p))^(1/p-1) * x_pde_hat_i{i}.^(p-1);
    end

    %% FILTERING AND MODIFICATION OF SENSITIVITIES
    dx_ind = cell(CellCount,1);
    for i = 1:CellCount
        dx_ind{i} = beta * (1-tanh(beta*(xTilde{i}-eta)).*tanh(beta*(xTilde{i}-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    end

    if FilterDom == 1
        dx = beta * (1-tanh(beta*(xTilde_cat-eta)).*tanh(beta*(xTilde_cat-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
        dc_i(:) = H*(dc_i(:).*dx(:)./Hs);
        %dv(:) = H*(dv(:).*dx(:)./Hs);
        dconstraint_temp=horzcat(dconstraints{1:CellCount});
        dconstraint_temp(:)=H*(dconstraint_temp(:).*dx(:)./Hs);
        for i=1:CellCount
            dconstraints{i}=dconstraint_temp(1:nely,(((i-1)*nelx)+1):i*nelx);
        end
    else
        for i = 1:CellCount
            dc_i_temp = dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx));
            dc_i_temp(:) = H_ind*(dc_i_temp(:).*dx_ind{i}(:)./Hs_ind);
            dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_i_temp;

            dconstraint_temp = dconstraints{i};
            dconstraint_temp(:) = H_ind*(dconstraint_temp(:).*dx_ind{i}(:)./Hs_ind);
            dconstraints{i} = dconstraint_temp;
        end
    end

    %% MMA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES

    m = (GLOBAL + LOCAL)*CellCount;
    n = CellCount*nelx*nely;

    move = 0.2;
    fval = zeros((GLOBAL+LOCAL)*CellCount,1);
    dfdx = zeros((GLOBAL+LOCAL)*CellCount,nelx*nely);

    c = c_i;
    f0val = c;
    Q=Q_i;
    df0dx = reshape(dc_i,[CellCount*nelx*nely,1]);
    %df0dx2 = 0*df0dx;

    M = 0;
    if GLOBAL == 1
        for i = 1:CellCount
            fval(M*CellCount+i) = 1000*constraints(i,1);
            dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))) = 1000*reshape(dconstraints{i,1},[1,nelx*nely]);
        end
        M = M + 1;
    end

    if LOCAL == 1
        for i = 1:CellCount
            fval(M*CellCount+i) = ((sum(x_pde_hat_i{i}.^p))^(1/p)- vol_max_pNorm(i));
            dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))) = TF'*(LF'\(LF\(TF*dfdx_pde_i{i}(:))));
            tmp = reshape(dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))),[nely,nelx]);
            dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))) = reshape(H_ind*(tmp(:).*dx_ind{i}(:)./Hs_ind),[1,nelx*nely]);
        end
        M = M + 1;
    end

    iter = loop;
    xval = reshape(x_cat,[CellCount*nelx*nely,1]);
    xmin = max(0.0,xval-move);
    xmax = min(1,xval+move);

    %dfdx2 = 0*dfdx;
    a0 = 1;
    a = zeros(m,1);
    c_ = ones(m,1)*1000;
    d = zeros(m,1);
    mdof = 1:m;

    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,fval(1:m),dfdx(mdof,:),low,upp,a0,a,c_,d);

    xnew_cat = reshape(xmma,[nely,CellCount*nelx]);
    xold2 = xold1;
    xold1 = xval;

    if FilterDom == 1
        xTilde_cat(:) = (H*xnew_cat(:))./Hs;
    else
        for i = 1:CellCount
            xTilde_temp = xTilde_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
            xnew_temp = xnew_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
            xTilde_temp(:) = (H_ind*xnew_temp(:))./Hs_ind;
            xTilde_cat(1:nely,(((i-1)*nelx)+1):(i*nelx)) = xTilde_temp;
        end
    end

    xPhys_i_cat = (tanh(beta*eta) + tanh(beta*(xTilde_cat-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    change = max(abs(xnew_cat(:)-x_cat(:)));
    x_cat = xnew_cat;

    for i = 1:CellCount
        x{i} = x_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
        xTilde{i} = xTilde_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
        xPhys_i{i} = xPhys_i_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));

        x{i} = ApplySymmetry(x{i},sym);
        xTilde{i} = ApplySymmetry(xTilde{i},sym);
        xPhys_i{i} = ApplySymmetry(xPhys_i{i},sym);
    end

    for i=1:CellCount
        xPhys_i_cat(1:nely,(((i-1)*nelx)+1):(i*nelx)) = xPhys_i{i};
    end

    i = CellCount + 1;
    l = 1;
    while i <= CaseCount
        j = 1;
        k = j + l;
        while k <= CellCount
            x{i} = horzcat(x{j},x{k});
            xTilde{i} = horzcat(xTilde{j},xTilde{k});
            xPhys_i{i} = horzcat(xPhys_i{j},xPhys_i{k});
            i  = i + 1;
            j = j + 1;
            k = j + l;
        end
        l = l + 1;
    end

    %% PRINT ITERATION PROCEDURES
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.4f\n',loop,c, ...
    mean(xPhys_i_cat(:)),change);
    fprintf(' Current Vol.:%7.3f Current Heat damp:%7.3f\n',current_vol,current_heatdamps);
    fprintf(' Cell target heat tensor constraints:\n');
    for i=1:CellCount
        fprintf(' Constraint %i :%7.3e\n',i,constraints(i));
    end

    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if beta < 256 && (loopbeta >= BetaUpdFreq || change <= 0.0001)
        beta = BetaUpdFactor * beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end

    %% PLOT DENSITIES
    figure(1);
    colormap(gray); imagesc(1-xPhys_i_cat); caxis([0 1]); axis equal; axis off; drawnow;
    plotadd=strcat(num2str(loop,'%04d'),'.png');
    imwrite(1-xPhys_i_cat,strcat('objective_history_demo\mutual,w=1\iteration_fig\',plotadd));
    xPhys_cat = xPhys_i_cat;

    %% PLOT OBJECTIVE FUNCTION
    if Influence == 1
        figure(2);
        clf;
        hold on
        for i = 1:CellCount
            plot(c_hist(1:loop,i),'b')
        end
    else
        figure(2);
        clf;
        subplot(1,4,[1,2])
        hold on
        for i = 1:CellCount
            plot(c_hist(1:loop,i),'b')
        end
        hold off
        title('Vols');
        
        subplot(1,4,[3,4])
        hold on
        for i = CellCount+1:CaseCount
            plot(c_hist(1:loop,i),'g')
        end
        hold off
        title('Heat Damps');
    end
    %autoArrangeFigures(2,2,1);
end
toc
end

%% FUNCTION: PREPARE FINITE ELEMENT ANALYSIS
function [nodenrs,edofMat,iK,jK] = PrepareFEA(nelx,nely)
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec=reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat=repmat(edofVec,1,4)+repmat([0,nely+[1,0],-1],nelx*nely,1);%指示每个有限元上分布的自由度下标，从左下角开始逆时针旋转。
iK=reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);%为了稀疏刚度矩阵准备的稀疏下标（行下标）
jK=reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);%为了稀疏刚度矩阵准备的稀疏下标（列下标）
end

%% FUNCTION: APPLY PERIODIC BOUNDARY CONDITIONS
function [T,Tfixed,d1,d2,d3,d4,wfixed] = ApplyPeriodicBC(nelx,nely,nodenrs)
e0=eye(2);%二阶单位矩阵，用途：每一列表示一种单元上的测试热流。分别是[1;0]和[0;1]
Tfixed=zeros(4,2);%对应原文中的T1，指代预设的角点温度。本文中选取矩形设计域的四个顶点上的四个自由度，注意有二种施加测试热流的情况需要考虑
T=zeros((nelx+1)*(nely+1),2);%所有自由度的温度的几何，注意有两种施加热流的场景，依次是：单位横向热流，单位纵向热流
alldofs=(1:(nelx+1)*(nely+1));%所有自由度的下标集合
n1=[nodenrs(end,[1,end]),nodenrs(1,[end,1])];%对应四个角点的下标，它们的温度是预设的，从左下角逆时针排列
d1=reshape(n1,1,4);%对应U1的自由度集合，即四个顶点上四个预设的自由度，从左下角逆时针排列
n3=[nodenrs(2:end-1,1)',nodenrs(end,2:end-1)];%左侧边界和下侧边界的格点下标集合，从左上角格点开始逆时针排列
d3=reshape(n3,1,nelx+nely-2);%对应U3的自由度集合，在U3上的自由度需要满足周期性边界条件，即和U4匹配
n4=[nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];%右侧边界和上侧边界的格点下标集合，需要和上文的D3一一对应
d4=reshape(n4,1,nelx+nely-2);%对应U4的自由度集合，注意需要和U3一一对应以方便设置边界条件
d2=setdiff(alldofs,[d1,d3,d4]);%所有无约束的内部自由度
%这下面的ufixed到底是怎么设定的？怎么和U3与U4产生的联系？现在判断，应该是和边界约束条件有关系，且默认左下角的两个自由度上的位移恒定为0
%边界对应关系是从左到右，从下到上的
for j=1:2
    Tfixed(2,j)=[e0(1,j),e0(2,j)]*[nelx;0];%右下角两个自由度的位移，同时也充当x向边界对应的差值关系
    Tfixed(4,j)=[e0(1,j),e0(2,j)]*[0;nely];%左上角两个自由度的位移，同时也中档y向边界对应的差值关系
    Tfixed(3,j)=Tfixed(2,j)+Tfixed(4,j);%右上角两个自由度的位移，同时考虑x和y向的差值关系
end
wfixed=[repmat(Tfixed(2,:),nely-1,1);repmat(Tfixed(4,:),nelx-1,1)];%对应原文中的W，用于控制边界条件
end

%% FUNCTION: INITIALIZE DESIGN VARIABLES
function [x, xTilde, xPhys_i, xold1, xold2] = InitializeIteration(nelx,nely,beta,eta,volfrac,IC,x_imp,imp)
if imp == 0
    %x = repmat(volfrac,nely,nelx);
    x = repmat(0.5,nely,nelx);
    for i = 1:nelx
        for j = 1:nely
            if IC == 1
                if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) < min(nelx,nely)/6
                    %x(j,i) = volfrac/2;
                    x(j,i) = 0;
                end
            elseif IC == 2
                if sqrt((i-nelx/4-0.5)^2+(j-nely/4-0.5)^2) < min(nelx,nely)/20
                    %x(j,i) = volfrac/2;
                    x(j,i) = 0;
                end
            elseif IC == 3
                if sqrt((i-nelx/4-0.5)^3+(j-nely/8-0.5)^2) < min(nelx/2,nely)/20
                    %x(j,i) = volfrac/2;
                    x(j,i) = 0;
                end
            elseif IC == 4
                if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) < min(nelx,nely)/20
                    %x(j,i) = volfrac/2;
                    x(j,i) = 0;
                end
                if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) < min(nelx/2,nely)/10 && sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) >= min(nelx/2,nely)/20
                    %x(j,i) = volfrac/2;
                    x(j,i) = 1;
                end
            elseif IC == 5
                if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) < min(nelx,nely)/8
                    %x(j,i) = volfrac/2;
                    x(j,i) = 0;
                end
                if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) < min(nelx/2,nely)/6 && sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) >= min(nelx/2,nely)/8
                    %x(j,i) = volfrac/2;
                    x(j,i) = 1;
                end
            end
        end
    end
else
    x = x_imp;
end

xTilde = x;
xPhys_i = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));

xold1 = reshape(x,[nely*nelx,1]);
xold2 = reshape(x,[nely*nelx,1]);
end

%% FUNCTION: FINITE ELEMENT ANALYSIS
function T = FEAnalysis(nelx,nely,KE,Kmin,xPhys,penal,K0,iK,jK,d1,d2,d3,d4,Tfixed,wfixed)
sK = reshape(KE(:)*(Kmin+xPhys(:)'.^penal*(K0-Kmin)),16*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
Kr = [K(d2,d2), K(d2,d3)+K(d2,d4); K(d3,d2)+K(d4,d2), K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
T(d1,:) = Tfixed;
T([d2,d3],:) = Kr\(-[K(d2,d1); K(d3,d1)+K(d4,d1)]*Tfixed-[K(d2,d4); ...
    K(d3,d4)+K(d4,d4)]*wfixed);
T(d4,:) = T(d3,:) + wfixed;
end

%% FUNCTION: HOMOGENIZATION
function [Q, dQ] = Homogenization(nelx,nely,T,edofMat,Kmin,K0,xPhys,penal,KE)
qe = cell(2,2);
Q = zeros(2);
dQ = cell(2,2);
for i=1:2
    for j=1:2
        U1=T(:,i);U2=T(:,j);
        qe{i,j}=reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx)/(nelx*nely);
        Q(i,j)=sum(sum((Kmin+xPhys.^penal*(K0-Kmin)).*qe{i,j}));
        dQ{i,j}=penal*(K0-Kmin)*xPhys.^(penal-1).*qe{i,j};
    end
end

end

%% FUNCTION: OBJECTIVE. CONDITION: 1,INDIVIDUAL 2,COMPOUND STRUCTURES
function [c,dc] = Objective(Q,dQ,OBJECTIVE,CONDITION)
if OBJECTIVE == 1
    c = 2/(Q(1,1)+Q(2,2));
    dc = -2*(dQ{1,1}+dQ{2,2})/((Q(1,1)+Q(2,2))^2);
elseif OBJECTIVE == 2
    if CONDITION == 1
        c = 1/Q(1,1);
        dc = -dQ{1,1}/(Q(1,1)^2);
    elseif CONDITION == 2
        c = 2/Q(1,1);
        dc = -2*dQ{1,1}/(Q(1,1)^2);
    end
elseif OBJECTIVE == 3
    if CONDITION == 1
        c = 1/Q(2,2);
        dc = -dQ{2,2}/(Q(2,2)^2);
    elseif CONDITION == 2
        c = 1/(2*Q(2,2));
        dc = -0.5*dQ{2,2}/(Q(2,2)^2);
    end
elseif OBJECTIVE == 4
    if CONDITION == 1
        c = 1/Q(1,1)+1/Q(2,2);
        dc = -dQ{1,1}/(Q(1,1)^2)-dQ{2,2}/(Q(2,2)^2);
    elseif CONDITION == 2
        c = 2/Q(1,1)+1/(2*Q(2,2));
        dc = -2*dQ{1,1}/(Q(1,1)^2)-0.5*dQ{2,2}/(Q(2,2)^2);
    end
end
end

%% FUNCTION: CONSTRAINTS. TO CONTROL HEAT TENSOR TARGETS
function [cons,dcons] = Constraint(Q,dQ,kalpa,epsilon)
%     cons=(Q(1,1)-kalpa(1)).^2/kalpa(1)+(Q(2,2)-kalpa(2)).^2/kalpa(2)-epsilon;
%     dcons=2*(Q(1,1)-kalpa(1))*dQ{1,1}/kalpa(1)+2*(Q(2,2)-kalpa(2))*dQ{2,2}/kalpa(2);
    cons=(Q(1,1)-kalpa(1)).^2+(Q(2,2)-kalpa(2)).^2-epsilon;
    dcons=2*(Q(1,1)-kalpa(1))*dQ{1,1}+2*(Q(2,2)-kalpa(2))*dQ{2,2};
end

%% APPLY SYMMETRY
function x = ApplySymmetry(x,sym)
nelx = size(x,2);

% 4-fold symmetry across diagonals
if sym == 1
    x = (x+x')/2;
    x = rot90(x);
    x = (x+x')/2;
    x = rot90(x);

% 4-fold symmetry across normals
elseif sym == 2
    for row = 1:nelx
        for col = 1:nelx/2
            x(row,col) = (x(row,col) + x(row,nelx-col+1))/2;
            x(row,nelx+1-col) = x(row,col);
        end
    end

    for col = 1:nelx
        for row = 1:nelx/2
            x(row,col) = (x(row,col) + x(nelx-row+1,col))/2;
            x(nelx+1-row,col) = x(row,col);
        end
    end

% 4-fold rotational symmetry
elseif sym == 3
      A = x;
      A90 = rot90(A);
      A180 = rot90(A,2);
      A270 = rot90(A,3);
      E = zeros(nelx);
      temp = 0.25*(A(1:nelx/2,1:nelx/2)+A90(1:nelx/2,1:nelx/2)+A180(1:nelx/2,1:nelx/2)+A270(1:nelx/2,1:nelx/2));
      for ii = 1:3
          E(1:nelx/2,1:nelx/2) = temp;
          E = rot90(E);
      end
      E(1:nelx/2,1:nelx/2) = temp;
      x = E;

% Cubic symmetry
elseif sym == 4
      A = x;
      A90 = rot90(A);
      A180 = rot90(A,2);
      A270 = rot90(A,3);
      E = zeros(nelx);
      temp = 0.25*(A(1:nelx/2,1:nelx/2)+A90(1:nelx/2,1:nelx/2)+A180(1:nelx/2,1:nelx/2)+A270(1:nelx/2,1:nelx/2));
      for ii = 1:3
          E(1:nelx/2,1:nelx/2) = temp;
          E = rot90(E);
      end
      E(1:nelx/2,1:nelx/2) = temp;
      x = E;

    for row = 1:nelx
        for col = 1:nelx/2
            x(row,col) = (x(row,col) + x(row,nelx-col+1))/2;
            x(row,nelx+1-col) = x(row,col);
        end
    end

    for col = 1:nelx
        for row = 1:nelx/2
            x(row,col) = (x(row,col) + x(nelx-row+1,col))/2;
            x(nelx+1-row,col) = x(row,col);
        end
    end
end
end
