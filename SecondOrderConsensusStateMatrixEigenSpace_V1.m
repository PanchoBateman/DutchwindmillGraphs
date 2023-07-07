%% This programm provides the eigenvalues and the eigenvector for a second order consensus dynamic matrix
% François Bateman 27/03/2022 
% Ecole de l'air et de l'Espace 
%
% Ltilde=[ 0 -I   ]
%        [ L xi*L ]
%
% L is the Laplacian matrix
% J is the Laplacian matrix in the Jordan Canonical Form
% V is the modal transformation matrix inv(V)*L*V=J
% xi is the damping factor
% I is the identity matrix
% Jtilde : Jordan Canonical Form of L
% Phi : modal transformation matrix  such that inv(Phi)*Ltilde*Phi=Jtilde
% 
% Graphs studied : Directed DutchWindMill_eta_kappa graphs
% eta : number of cycles
% kappa : number of edges for each cycle
% 
%%
% Use a precomputed graph 



clear all 
clc
display('Open a file containing the Adjacency, Laplacian, Jordan Canonical Form and the associated Transformation Matrix for the eta kappa DutchwindMill graph') 
file=uigetfile
load(file)





% Display the graph 
 G=digraph(A) 
 figure(1)
 plot(G)
 title(['Dutchwindmill',num2str(eta),num2str(kappa)]);
% Eigenvalues of L
 Lambda=diag(J);
 n=length(J);
% Damping factor
 xi=1;
% Second order dynamic consensus matrix
 Ltilde=[zeros(n,n),-eye(n);L xi*L];
 
 
 
% Augmented Laplacian matrix Ltild eigenvalues
 
 for j=1:n;
     mu(j)=(xi*Lambda(j)+sqrt((xi*Lambda(j))^2-4*Lambda(j)))/2;
     mu(j+length(L))=(xi*Lambda(j)-sqrt((xi*Lambda(j))^2-4*Lambda(j)))/2;      
 end
 mu=transpose(mu);
 
 V=[V,V];
 
 % Jordan Blocks indice
 colsingle=1:n; 
  for j=1:n-1;
    updiaglambda(j+1)=J(j,j+1); 
    updiaglambda=round(updiaglambda,10); % modifiée avec round pour utiliser ComputeDWGraphJCF (different de 1 et complexe sur la supérieure de la diagonale)
    coljordan= [find(updiaglambda==1)];  % index one corresponds with the zero eigenvalue
    
  end
 
 colsingle(coljordan)=[];
 
 Phi_=zeros(2*n,2*n);
 
 %% Ltild eigenvectors calculation

% Eigenvectors associated with the two zero eigenvalues (columns 1 & 1+n)
Phi_(:,1)=[1/sqrt(n)*ones(n,1) ; zeros(n,1)];
 %Phi(:,1+n)=[zeros(n,1) ; -1/sqrt(n)*ones(n,1)];
Phi_(:,1+n)=[1/sqrt(n)*ones(n,1) ; -1/sqrt(n)*ones(n,1)];
 
% Eigenvectors associated with the single nonzero eigenvalues
  
 for j=2:length(colsingle),
   Phi_(:,colsingle(j))=  [V(:,colsingle(j)); -mu(colsingle(j))*V(:,colsingle(j))]
   Phi_(:,colsingle(j)+n)=  [V(:,colsingle(j)+n); -mu(colsingle(j)+n)*V(:,colsingle(j)+n)]
 end
 
% Eigenvectors associated with eigenvalues whose algebraic mul >geometric
% The eigenvalues are stored in an array, each column
mu1=reshape(mu(coljordan),sizejordan,nb_jordan);
mu2=reshape(mu(coljordan+n),sizejordan,nb_jordan);
 
% Redimmensionnement du tableau des indices colonnes de la matrice de
% passage V associés aux blocs de Jordan
% row : dim d'un bloc de Jordan -1
% column : nombre de blocs de Jordan
blocjordan=reshape(coljordan,sizejordan,nb_jordan);
 
 
  
%% Ltild eigenvectors calculation


% Computation of the transformation matrix Phi of Ltilde

alpha=zeros(sizejordan,sizejordan+1,nb_jordan); 
beta=zeros(sizejordan,sizejordan+1,nb_jordan); 
%Phi_=[];



%% 

for j=1:nb_jordan, 
     i=1;
     d=(mu1(i,j)-mu2(i,j));
     q=(mu2(i,j)/mu1(i,j));
     xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
     lambda=mu1(i,j)*mu2(i,j);
     alpha(i,i+1,j)= (d*q)^i;
       
    
     for i=2:sizejordan,   
       d=(mu1(i,j)-mu2(i,j));
       q=(mu2(i,j)/mu1(i,j));
       xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
       lambda=mu1(i,j)*mu2(i,j);
             
       %beta(i+1,i)=alpha(i,i);
       alpha(i,i+1,j)= (d*q)^i;
       for k=2:i,       
            if i==2 & k==2,
                  beta(i-1,k-1,j)=1;
                  beta(i,k,j)=d*q;    
            end
           alpha(i,k,j)=q*(d*alpha(i-1,k-1,j)+beta(i-1,k-1,j)-xi*alpha(i-1,k,j)) 
           if i==2 & k==2,
                  beta(i-1,k-1,j)=1;
           else
                beta(i,k,j)= alpha(i-1,k,j);
           end
       end
       
      
     end
end

alpha(:,1,:)=[] ; % suppress the first null columns of arrays alpha
for j=1:nb_jordan,
    for i=1:sizejordan
        
        for k=1:sizejordan,
            blocjordan(i,j)
            Phi_(1:n,blocjordan(i,j))=alpha(i,k,j)*V(:,blocjordan(k,j))+Phi_(1:n,blocjordan(i,j));
            Phi_(n+1:2*n,blocjordan(i,j))=-Phi_(1:n,blocjordan(i,j)-1)-mu1(i,j)*Phi_(1:n,blocjordan(i,j));
           
        end
       
    end
end

%%

clear alpha beta
for j=1:nb_jordan, 
     i=1;
     d=(mu2(i,j)-mu1(i,j));
     q=(mu1(i,j)/mu2(i,j));
     xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
     lambda=mu1(i,j)*mu2(i,j);
     alpha(i,i+1,j)= (d*q)^i;
     
    
    
    
     for i=2:sizejordan,   
       d=(mu2(i,j)-mu1(i,j));
       q=(mu1(i,j)/mu2(i,j));
       xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
       lambda=mu1(i,j)*mu2(i,j);
             
       %beta(i+1,i)=alpha(i,i);
       alpha(i,i+1,j)= (d*q)^i;
       for k=2:i,       
            if i==2 & k==2,
                  beta(i-1,k-1,j)=1;
                  beta(i,k,j)=d*q;    
            end
           alpha(i,k,j)=q*(d*alpha(i-1,k-1,j)+beta(i-1,k-1,j)-xi*alpha(i-1,k,j))
           if i==2 & k==2,
                  beta(i-1,k-1,j)=1;
           else
                beta(i,k,j)= alpha(i-1,k,j);
           end
       end
       
      
     end
end
%Phi_=zeros(n,n);
alpha(:,1,:)=[] ; % suppress the first null columns of arrays alpha
for j=1:nb_jordan,
    for i=1:sizejordan
        
        for k=1:sizejordan,
            Phi_(1:n,blocjordan(i,j)+n)=alpha(i,k,j)*V(:,blocjordan(k,j)+n)+Phi_(1:n,blocjordan(i,j)+n);      
            Phi_(n+1:2*n,blocjordan(i,j)+n)=-Phi_(1:n,blocjordan(i,j)+n-1)-mu2(i,j)*Phi_(1:n,blocjordan(i,j)+n);
          
        end
        %
    end
end



%% Check that Phi is the transformation matrix for the JCF Jtilde of Ltilde
Jtilde= inv(Phi_)*Ltilde*Phi_

figure(2)
spy(round(Jtilde,10))
figure(3)
spy(round(Phi_,10))
 


