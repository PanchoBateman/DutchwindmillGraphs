function [A,L,J,P,eta, kappa,nb_jordan,sizejordan,coljordan]=ComputeDWMGraphsJCF
% This program computes the eigenvalues and the eigenvectors of the
% Laplacian matrix L for a directed Dutch windmill (DWM) graph.
% The DWM graph is compound with eta cycles and each cycles has kappa
% vertices + a common vertex. 
% To display the Jordan chains : uncomment lines 130-132
%% Syntax : [A,L,J,P,eta, kappa,nb_jordan,sizejordan,coljordan]=ComputeDWMGraphsJCF_V1
% A : adjacency matrix
% L : Laplacian matrix
% J the Jordan Canonical Form of L
% P the transformation matrix 
% eta : number of cycles
% kappa : each cycle has kappa vertices + the common vertex
% nb_jordan : nuber of Jordan blocks
% sizejordan : size of each Jordan block 
% coljordan : vector containing the indices of the non diagoan "1" associated with the Jordan blocks 

 
% François BATEMAN 09/06/2023
% Centre de Recherche de l'Ecole de l'Air (et de l'Espace)
% francois.bateman@ecole-air.fr

%%

clear all
close all
clc
digits(32)

disp('---------------------------------------------------------------------------------------------------------------------------------------')      
disp('Compute the Adjacency, Laplacian, Jordan Canonical Form and the associated Transformation Matrix for the eta kappa DutchwindMill graph '), 
disp('---------------------------------------------------------------------------------------------------------------------------------------')    
disp('---------------------------------------------------------------------------------------------------------------------------------------')    

% Define the graph
if nargin<1 | nargin>3 
    eta=input('Number of cycles eta : ');
    kappa=input('Number of nodes for each cycle kappa : ');
end

     
n=eta*kappa+1;
i=[1:eta*kappa];
s=reshape(i,kappa,eta);
s(kappa+1,:)=zeros(1,eta);
s=reshape(s,1,kappa*eta+eta);
t=wshift('1D',s,-1) ;
s=s+1;
t=t+1;
G=digraph(t,s);
A=full(adjacency(G));
D=diag(indegree(G));
L=D-A;            


% Display the graph 
names={'0'};
for k=1:eta*kappa+1,
    names{k}=num2str(k-1);
end 
G=digraph(A',names) ;
%clear names
figure(1)
plot(G)
clear k name names
title(['Dutchwindmill',' \eta= ',num2str(eta),' \kappa= ',num2str(kappa)]);


%% Special Case eta = 4 and kappa = 2
if eta==4 & kappa==2,
[P,J]=jordan(L);
lambda=diag(J);
nb_jordan=4;
sizejordan=2;
V=P;

else
 
%% Initialization of P, V, V0 W 
P=[]; V=[]; V0=[]; W=[];

for l=1:kappa,
            
    % Computes the Jordan's chains (L-1I)^l
    e1=zeros(kappa,1); e1(1)=1;
    ek=zeros(kappa,1);ek(kappa)=1;
    L1=(L-eye(n))^l;


    %% Builds the (L-I]^l matrix from equation (11)

    M=zeros(kappa,kappa);
    N=M;
    for i=1:kappa-1,
        N(i,i+1)=1;
    end


    for j=1: l-1,
        for i=0:j-1,
            M=M+(-1)^(j-1)*(eta-1)^(l-j-1)*N^(j-i-1)*ek*e1.'*N^i;
        end
    end

    R=zeros(kappa,1)';

    for j=1:l,
        R=R+ (-1)^j*(eta-1)^(l-j)*e1.'*N^(j-1);
    end


    S=zeros(kappa,1)';
    for j=1:l,
        e=zeros(kappa,1);
        e(kappa-j+1)=1;
        S=(-1)^j*(eta-1)^(l-j)*transpose(e)+S;
    end
    S=S';

    L_1=zeros(n,n);
    L_1(1,1)=(eta-1)^l;

    for i=2:kappa:(eta-1)*kappa+2,  
        for j=2:kappa:(eta-1)*kappa+2;
            L_1(i:i+kappa-1,j:j+kappa-1)=M ;    
        end
    end    

    for j=2:kappa:(eta-1)*kappa+2,  
        L_1(1,j:j+kappa-1)=R;
        L_1(j:j+kappa-1,1)=S;
        L_1(j:j+kappa-1,j:j+kappa-1)=M+(-1)^l*(N^l);  
    end
    L_1(:,1+l)  ;  

% Uncomment these lines to display the Jordan's chains    
% display(['(L-In)^',num2str(l),'= ']); 
% L_1
% pause(1);


    
%% Computation of the eigenvectors for the 1 eigenvalue

   
    p=zeros(eta*kappa+1,1);
   
        for j=l+1:kappa:eta*kappa+1,          
            p(j)=1;                     % build the eigenstructure
        end
    nullp=null(transpose(p),'r');       %%%%%% Là est peut être le pb de conditionnement
    %nullp=null(transpose(p));           % transpose of  the nullspace of p and provides the desired eigenvectors + undesired vectors (inapropriate structure)
    p=(p(1:end,1).*nullp(1:end,:));     % suppress the vectors of the nullspace which structure is inapropriate
    
    for i=1:(eta)*kappa,               % we keep eta-1 eigenvectors among eta   
    %for i=1:(eta-1)*kappa,               
        if ~isempty(find(p(:,i)~=0)),
         
        if mod(i,2)==0,      
           P=[P,p(:,i)]  ;                % P contains the (eta-1)*kappa eigenvectors associated with the eigenvalues {1}
        else                 
           P=[P,-p(:,i)]  ;   %ajout pour avoir +1 au lieu de -1 dans la matrice de jordan 
         
        end                  
     end
    
   end
   

    
  
    
end

%% Reordering of the colums of the transformation matrix P    
   indice=1:(eta-1)*kappa;                    
   indice=reshape(indice,eta-1,kappa );
   P=P(:,reshape(indice',(eta-1)*kappa,1));


%% Computation of the others eigenvalues and eigenvectors p0=1

% roots of the polynomial (eta-x)(1-x)^kappa-eta
% C=[];
% for k=0:kappa,
%     C=[C,nchoosek(kappa,k)*(-1)^k]; % coefficients of (1-x)^eta
% end
%      C=conv([-1 eta],[C]);
%      C(:,kappa+2)=0;
%      lambda=roots(C);               % non unit eigenvalues of L
%      lambda=sort(lambda);           % sort the eigenvalues (Ajout tardif pas utile)
% V0=ones(1,kappa+1);
% for k=1:kappa+1,                    % compute the eigenvectors associated with the previous eigenvalues
% %V((k-1)*kappa*eta+1)=1    
%      for j=1:eta,   
%          v(1)=(eta-lambda(k))/eta;
%          for m=2:kappa,             
%              v(m)=(1-lambda(k))*v(m-1);
%          end
%      v=transpose(v); 
%      V=[V;v];
%      v=[];
%      end
%     
% end
% V=reshape(V,eta*kappa,kappa+1);
% 
% V=[V0;V];

%% Computation of the others eigenvalues and eigenvectors pk=1 and Vo=1/n*(1 1 1...1)'

%roots of the polynomial (eta-x)(1-x)^kappa-eta

C=[];
for k=0:kappa,
    C=[C,nchoosek(kappa,k)*(-1)^k]; % coefficients of (1-x)^eta
end
     C=conv([-1 eta],[C]);
     C(:,kappa+2)=0;
     lambda=roots(C);               % non unit eigenvalues of L
     lambda=sort(lambda);           % sort the eigenvalues (Ajout tardif pas utile)
V0=ones(eta*kappa+1,1);
for k=2:kappa+1,                    % compute the eigenvectors associated with the previous eigenvalues
   
     for j=1:eta,   
         v(kappa)=1;
         for m=kappa:-1:2,             
             v(m-1)=v(m)/(1-lambda(k));
         end
         v=transpose(v); 
         V=[V;v];
         v=[];
         w=1-lambda(k);
         
     end
  W=[W,w];  
end
V=reshape(V,eta*kappa,kappa);
V=[W;V]; % original

V=[V0,V];


%% Computes the transformation matrix P

%P=[V P];
P=[V(:,1),P,V(:,2:end)]; % original

%% Computes the JCF J

J=inv(P)*L*P;




%V=P;
nb_jordan=eta-1;
sizejordan=kappa-1;
% eta=eta;
% kappa=kappa;

n=length(J);
 colsingle=1:n; 
  for j=1:n-1;
    updiaglambda(j+1)=J(j,j+1); 
    updiaglambda=round(updiaglambda,10); % modifiée avec round pour utiliser ComputeDWGraphJCF (different de 1 et complexe sur la supérieure de la diagonale)
    coljordan= [find(updiaglambda==1)];  % index one corresponds with the zero eigenvalue
    
  end
  
% clean and display J


%  J=diag(diag(J));
%  J(1,1)=0;
% for k=1:length(coljordan)
%     J(coljordan(k)-1,coljordan(k))=1;
% end


end  %(case eta=4 and kappa=2)



%% Display the structures of matrices J and P


figure(2)
subplot(121)
%spy(round(abs(J),p))
spy(round(abs(J)))
title('Structure of the JCF of the J matrix')
subplot(122)
%spy(round(abs(P),q))
spy(round(abs(P)))
title('Structure of the transformation matrix P')
sgtitle(['Directed Dutchwindmill',' \eta= ',num2str(eta),' \kappa= ',num2str(kappa)] ) 

%% Backup 

filename=['Dutchwindmill',num2str(eta),num2str(kappa),'.mat'];
reply=input(['Save these results in the ',filename, ' file ? [y/n] '],'s');
save(filename)
end