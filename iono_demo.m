
load iono;


% rand_sequence=randperm(size(iono,1));
% temp_dataset=iono;
% iono=temp_dataset(rand_sequence, :);


P=iono(1:200,2:(size(iono,2)))';
T=iono(1:200,1)';
pos5=find(T==0);

T(pos5)=-1;

xapp=P';
yapp=T';

P1=iono(201:size(iono,1),2:(size(iono,2)))';
T1=iono(201:size(iono,1),1)';
pos6=find(T1==0);
T1(pos6)=-1;

% P=mapminmax(P);
% P1=mapminmax(P1);


xtest=P1';
ytest=T1';


dim=size(P,1);

 % C=a1;
 
c = 1000;
 
 
  kerneloption=rand(1,dim)*dim;



%-----------------------------------------------------
%   Learning and Learning Parameters

% c =a1(aa1);
% kerneloption=a2(aa2);


  % c =a1(aa1);
 
 
   kerneloption=rand(1,dim)*dim;




epsilon = .000001;
kernel='gaussian';
verbose = 0;
tic
[xsup,w,b,pos,pos1,pos2,pos3,ps,H,n3]=svmclassrand(xapp,yapp,c,epsilon,kernel,kerneloption,verbose);
t=toc;


%--------------Train data error-------------------------%
y= svmvalrand(xapp,xsup,w,b,kernel,kerneloption,n3);
%---------------Accuracy of trian data-----------------------%

   
   tic
    y1= svmvalrand(xtest,xsup,w,b,kernel,kerneloption,n3);
    t1=toc;
    

   n1=size(P1,2);
     err1=0;
   for i=1:n1
    if(T1(i)==1)
        if(y1(i)<=0)
            err1=err1+1;
        end
    else
        if(y1(i)>=0)
            err1=err1+1;
        end
    end
   end
   
   accu1=(n1-err1)/n1

