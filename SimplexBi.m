function [ Beff,xExtremEff ] = simplexBi( type,c1,c2,A,b )

if (type == 'max')
  c1=-c1;
  c2=-c2;
end

[m,n] = size(A) ; 
%ajout des variable d'écart 
A = [A eye(m)];
b = b(:);
c1 = c1(:)';
c2=c2(:)'; 
A = [A b];
d = [c1 zeros(1,m+1);c2 zeros(1,m+1)];
beff=[];  % ensemble des bases efficientes 
xeff=[]; %ensembles des points extremes efficients 
zdom=[] ; %ensembles des vecteurs dominants 
A= [A;d];

subs = n+1:m+n; %base réalisable 
 
 lambda= 1 ;

  clambda= c1*lambda + (1-lambda)*c2 ; %calcul de clambda 
  disp(sprintf('     ________________________________________________'))
   disp(sprintf('\n              Tableaux de Simplexe '))
   disp(sprintf('     ________________________________________________'))
   disp(sprintf('\n                  Tableau Initiale\n'))
   A
   disp(sprintf(' Cliquer sur une touche pour continuer...\n\n'))
   pause
  while 1==1 
     clambda= A(m+1,1:n+m)*lambda + (1-lambda)*A(m+2,1:n+m) ;
     disp('Clambda'); disp(clambda);
    
   while any(clambda<-eps)   % Si un élément de clambda est négatif, la base n'est pas efficiente

       
       mi=min(clambda) ; 
       col = find(clambda==mi);
       disp('this') ; disp(col) ;
       [row,ratio] = MRT(A(1:m,m+n+1),A(1:m,col));
       if ~isempty(row)
      A(row,:)= A(row,:)/A(row,col);
      subs(row) = col;
      for i = 1:m+2
         if i ~= row
            A(i,:)= A(i,:)-A(i,col)*A(row,:);
         end
      end
       end
              disp(A);
                     clambda= A(m+1,1:n+m)*lambda + (1-lambda)*A(m+2,1:n+m) ;
  end
         % Clambda est positif donc la base est efficiente 
         z1 = A(m+1,m+n+1);
      z2 =A(m+2,m+n+1);
      zdom=[zdom ; z1 z2] ; 
      x = zeros(1,m+n);
      x(subs) = A(1:m,m+n+1);
      x = x(1:n);
        beff=[beff; subs]; 
        xeff=[xeff;x] ;
        I=find(A(m+2,1:m+n)<0); 
       if ~isempty(I)
        [lambda,ind]=max(-A(m+2,I)/(A(m+1,I)-A(m+2,I))); 
         bb=I(ind); 
         
         [row,ratio] = MRT(A(1:m,m+n+1),A(1:m,bb));
        
       if ~isempty(row)
      A(row,:)= A(row,:)/A(row,bb);
      subs(row) = bb;
      for i = 1:m+2
         if i ~= row
            A(i,:)= A(i,:)-A(i,bb)*A(row,:);
         end
      end
       end
       disp(A);
       else
           disp('Terminé...');
           break;
       end
        pause;
 
  end
  disp('Les solutions extremes efficientes sont : ');
  disp(xeff); 
  disp('Les bases efficientes sont : ');
  disp(beff);
  disp('Les points des vecteurs critères non dominés sont : ') ;
  disp(zdom); 
    
       
       
 

 
end
