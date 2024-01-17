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
%  Caution : The current version of this program does not deal with the case where the damping factor xi is equal to 2. 
% The program functions for xi<xi_min but the consensus is not reach.
%%
% Use a precomputed graph 
clc
clear all 
close all
disp('________________________________________________________________________________________________________')
disp('| This program computes the Jordan Canonical Form of the second order consensus state space matrix      |')
disp('|                                                                                                       |');
disp('|                                                                                                       |');
disp('|       Ltilde=[ 0 -I   ]                                                                               |')
disp('|              [ L xi*L ]                                                                               |')
disp('|                                                                                                       |');
disp('|                                                                                                       |');
disp('| when the algebraic multiplicity of the Laplacian matrix is greater than their geometric multiplicity  |') 
disp('|                                                                                                       |')
disp('| Caution : for the current release of this programm, xi must be different of 2                         |')
disp('|                                                                                                       |');
disp('|_______________________________________________________________________________________________________|')
disp(' ');
flag=input('Load a precomputed graph (1) or define and compute a graph (2) : ');


switch flag
    case 1
    % Display a pre-computed Laplacian matrix 
      display('Open a file containing the Adjacency, Laplacian, Jordan Canonical Form and the associated Transformation Matrix for the eta kappa DutchwindMill graph '), 
      file=uigetfile;
      load(file)
    
     % Display the graph 
     % REMPLACER eta*kappa+1 PAR n
      names={'0'};
      for k=1:n%,eta*kappa+1,
          names{k}=num2str(k-1);
      end 
      G=digraph(A',names) ;
      figure(1)
      plot(G)
      clear k names
      %title(['Dutchwindmill',' \eta= ',num2str(eta),' \kappa= ',num2str(kappa)]);
    case 2
    % Given eta, kappa, this part of the program computes the Laplacian matrix 
    [A,L,J,P,eta, kappa,nb_jordan,sizejordan]=ComputeDWMGraphsJCF;
end





% Eigenvalues of L
 Lambda=diag(J);
 n=length(J);
% Damping factor
% Condition for 2nd order consensus on the damping factor - ref biblio [9]
xi_min= sqrt(max(((abs(imag(Lambda))).^2)./((abs(real(Lambda))).*abs(Lambda).^2)));
disp(['For this graph, to reach the consensus, the damping factor xi must be greater than ',num2str(xi_min)])
xi=input('Damping factor xi (default value 1) = ')
if isempty(xi) | xi==2,
    xi=1;
end
 
 
% Second order dynamic consensus matrix
 Ltilde=[zeros(n,n),-eye(n);L xi*L];
 
 
 
% Augmented Laplacian matrix Ltild eigenvalues
 
 for j=1:n,
     mu(j)=(xi*Lambda(j)+sqrt((xi*Lambda(j))^2-4*Lambda(j)))/2;
     mu(j+length(L))=(xi*Lambda(j)-sqrt((xi*Lambda(j))^2-4*Lambda(j)))/2;      
 end
 mu=transpose(mu);
 
 V=[P,P];
 
 % Jordan Blocks indice
 colsingle=1:n; 
  for j=1:n-1,
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
   Phi_(:,colsingle(j))=  [V(:,colsingle(j)); -mu(colsingle(j))*V(:,colsingle(j))];
   Phi_(:,colsingle(j)+n)=  [V(:,colsingle(j)+n); -mu(colsingle(j)+n)*V(:,colsingle(j)+n)];
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




%% 

for j=1:nb_jordan, 
     i=1;
     d=(mu1(i,j)-mu2(i,j));
     q=(mu2(i,j)/mu1(i,j));
     %xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
     lambda=mu1(i,j)*mu2(i,j);
     alpha(i,i+1,j)= (d*q)^i;
       
    
     for i=2:sizejordan,   
       d=(mu1(i,j)-mu2(i,j));
       q=(mu2(i,j)/mu1(i,j));
       %xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
       lambda=mu1(i,j)*mu2(i,j);
             
       %beta(i+1,i)=alpha(i,i);
       alpha(i,i+1,j)= (d*q)^i;
       for k=2:i,       
            if i==2 & k==2,
                  beta(i-1,k-1,j)=1;
                  beta(i,k,j)=d*q;    
            end
           alpha(i,k,j)=q*(d*alpha(i-1,k-1,j)+beta(i-1,k-1,j)-xi*alpha(i-1,k,j)) ;
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
    for i=1:sizejordan,
        
        for k=1:sizejordan,
            blocjordan(i,j);
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
     %xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
     lambda=mu1(i,j)*mu2(i,j);
     alpha(i,i+1,j)= (d*q)^i;
     
    
    
    
     for i=2:sizejordan,   
       d=(mu2(i,j)-mu1(i,j));
       q=(mu1(i,j)/mu2(i,j));
       %xi=((mu1(i,j)+mu2(i,j))/(mu1(i,j)*mu2(i,j)));
       lambda=mu1(i,j)*mu2(i,j);
             
       %beta(i+1,i)=alpha(i,i);
       alpha(i,i+1,j)= (d*q)^i;
       for k=2:i,       
            if i==2 & k==2,
                  beta(i-1,k-1,j)=1;
                  beta(i,k,j)=d*q;    
            end
           alpha(i,k,j)=q*(d*alpha(i-1,k-1,j)+beta(i-1,k-1,j)-xi*alpha(i-1,k,j));
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
    for i=1:sizejordan,
        
        for k=1:sizejordan,
            Phi_(1:n,blocjordan(i,j)+n)=alpha(i,k,j)*V(:,blocjordan(k,j)+n)+Phi_(1:n,blocjordan(i,j)+n);      
            Phi_(n+1:2*n,blocjordan(i,j)+n)=-Phi_(1:n,blocjordan(i,j)+n-1)-mu2(i,j)*Phi_(1:n,blocjordan(i,j)+n);
          
        end
        
    end
end



%% Check that Phi is the transformation matrix for the JCF Jtilde of Ltilde
Jtilde= inv(Phi_)*Ltilde*Phi_;


%% Display the structures of matrices J and P
% For large graphs, this part slows down the execution of the program

x=min(abs(Lambda(2:end)));  %find the smallest non zero eigenvalue magnitude
% number of digits after the decimal point  https://fr.mathworks.com/matlabcentral/profile/authors/2872967

p=0;
q=0;


while (floor(x*10^p)<100),%~=x*10^p)
p = p+1;
 
end


y=min(abs(P(abs(P)>0)));
while (floor(y*10^q)~=y*10^q)
q=q+1;
end



[row_0,col_0]=find(round(abs(P),q+3)==0);
for u=1:length(row_0)
    P(row_0(u),col_0(u))=0;
end

figure(3)

subplot(121)
spy(round(10*abs(Jtilde),p))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for the paper
% cmap=colormap('gray');
% size(cmap)
% cmap=flip(cmap,1);
% spyc(floor(10*abs(Jtilde)),cmap,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



title('Structure of the JCF $\tilde{J}$ of the state space matrix', 'Interpreter', 'LaTeX')

subplot(122)
spy(round(abs(Phi_),q))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for the paper
% cmap=colormap('copper');
% size(cmap)
% cmap=flip(cmap,1);
% spyc(floor(10*abs(Phi_)),cmap,0);
% set(gca, 'YTickLabel',{'','x_0','x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8','x_9','x_{10}','x_{11}','x_{12}','x_{13}','x_{14}','x_{15}','x_{16}''v_0','v_1','v_2','v_3','v_4','v_5','v_6','v_7','v_8','v_9','v_{10}','v_{11}','v_{12}','v_{13}','v_{14}','v_{15}'})
% set(gca, 'XTickLabel',{'','\mu_{01}','\mu_{11}','\mu_{21}','\mu_{31}','\mu_{41}','\mu_{51}','\mu_{61}','\mu_{71}','\mu_{71}','\mu_{71}','\mu_{71}','\mu_{71}','\mu_{71}','\mu_{71}','\mu_{71}','','\mu_{02}','\mu_{12}','\mu_{22}','\mu_{32}','\mu_{42}','\mu_{52}','\mu_{62}','\mu_{72}','\mu_{72}','\mu_{72}','\mu_{72}','\mu_{72}','\mu_{72}','\mu_{72}','\mu_{72}'})
% title('Structure of the transformation matrix $\Phi$','Interpreter', 'LaTeX')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sgtitle('$2^{nd}$ order consensus : JCF of the state space matrix ' ,'Interpreter', 'LaTeX')