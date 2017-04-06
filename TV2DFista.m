function [X_den,iter,fun_all]=TV2DFista(Xobs,lambda,l,u,pars)
%This function implements the FISTA method for TV-based denoising problems
%
% Based on the paper
% Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
% Total Variation Image Denoising and Deblurring Problems"
% -----------------------------------------------------------------------
% Copyright (2008): Amir Beck and Marc Teboulle
% Modified by Nick Antipa, UC Berkeley, 2017
%
% FISTA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% INPUT
% Xobs ..............................an observed noisy image.
% lambda ........................ parameter
% l ..................................... lower bound on the pixels' values
% u ..................................... upper bound on the pixels' values
% pars.................................parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.epsilon ..................... tolerance for relative error used in
%                                                       the stopping criteria (Default=1e-4)
% 
% OUTPUT
% X_den ........................... The solution of the problem 
%                                            min{||X-Xobs||^2+2*lambda*TV(X
%                                            ) : l <= X_{i,j} <=u} 
% iter .............................  Number of iterations required to get
%                                            an optimal solution (up to a tolerance)
% fun_all ......................   An array containing all the function
%                                             values obtained during the
%                                             iterations


%Define the Projection onto the box

% Assigning parameres according to pars and/or default values
flag=exist('pars','var');
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if (flag&&isfield(pars,'epsilon'))
    epsilon=pars.epsilon;
else
    epsilon=1e-4;
end

project = @(x)min(max(x,l),u);
Ltv = @(P1,P2)cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));
[m,n]=size(Xobs);
R1=zeros(m-1,n);
R2=zeros(m,n-1);
P1=zeros(m-1,n);
P2=zeros(m,n-1);
tkp1=1;
count=0;
i=0;

D=zeros(m,n);
fun_all=zeros(1,MAXITER);

while((i<MAXITER)&&(count<5))
    %%%%%%%%%
    % updating the iteration counter
    i=i+1;
    %%%%%%%%%
    % Storing the old value of the current solution
    Dold=D;
    %%%%%%%%%%
    %Computing the gradient of the objective function
    Pold1 = P1;
    Pold2 = P2;
    tk=tkp1;
    D=project(Xobs-lambda*Ltv(R1,R2));
    Q1 = -diff(D,1,1);
    Q2 = -diff(D,1,2);
    
    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    P1=R1+1/(8*lambda)*Q1;
    P2=R2+1/(8*lambda)*Q2;
    
    %%%%%%%%%%
    % Peforming the projection step

    A=[P1;zeros(1,n)].^2+[P2,zeros(m,1)].^2;
    A=sqrt(max(A,1));
    P1=P1./A(1:m-1,:);
    P2=P2./A(:,1:n-1);


    %%%%%%%%%%
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    R1=P1+(tk-1)/tkp1*(P1-Pold1);
    R2=P2+(tk-1)/tkp1*(P2-Pold2);
    
    re=norm(D-Dold,'fro')/norm(D,'fro');
    if (re<epsilon)
        count=count+1;
    else
        count=0;
    end
    
    C=Xobs-lambda*Ltv(P1,P2);
    PC=project(C);
    fval=-norm(C-PC,'fro')^2+norm(C,'fro')^2;
    fun_all(i) = fval;

end
fun_all = fun_all(1:i);
X_den=D;
iter=i;

