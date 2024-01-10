clear all;
tic;

 %%

filename = 'poly1_05'; %%%%Change file name when required...
fileread = sprintf('%s.dat',filename);
N = dlmread(fileread);
ScF=1;

N(:,1:4) = ScF*N(:,1:4);
[a0w,b0w] = max(N(:,3));
[a0m,b0m] = min(N(:,3));
Wall=[-0.0125,0.0125;-0.0125,0.0125;a0m-N(b0m,4),a0w+N(b0w,4)];

fprintf('Particle DATA is read\n\n');

Xi = 0.71;

Top_temp = 10;
Bot_temp = 20;

Delta_T = abs(Top_temp-Bot_temp);

%%
Tc = 25;
Tk = Tc+273.14;

% K_s=100;
%Li4OSi4 
K_s = (((7.317E-12)*(Tc.^(4.0)))+((-1.302E-8)*(Tc.^(3.0)))+(8.712E-6*(Tc.^(2.0)))+(-0.002876*Tc)+2.620);
 
%Li2TiO3
% poro = 0.1;
% kbeta = 1.06-(2.88E-4)*Tk;
% K_s = ((1-poro)/(1+kbeta*poro))*(4.77 - (5.11E-3)*Tk + (3.12E-6)*Tk*Tk);

%%
counter_file = 0;

dir_other = [2,3;1,3;1,2];
File_part = zeros(1,20);

%%
%Change the below for loop for specific direction
for direction_cond = 3:1:3

other_dir = dir_other(direction_cond,:);   


%Sorting through the columns based on direction
N = sortrows(N,direction_cond);

N_particles = size(N,1);
mean_r  = mean(N(:,4));

N(:,8) = 1:1:N_particles;

%%
%Creating the contact list
C_List = contact_list_va(N,0,direction_cond);

N_contacts = size(C_List,1);

fprintf('Contact list is generated for direction %d\n',direction_cond);

Cond = zeros(N_contacts,1);

N_overlap=0;

%Conduction Matrix

for i=1:1:N_contacts
 
    P_1 = C_List(i,1); 
    P_2 = C_List(i,2);   
    
    C_List(i,10) = N(P_1,direction_cond);
    C_List(i,11) = N(P_2,direction_cond);
    
    d = C_List(i,3);
    r_1 = N(P_1,4);
    r_2 = N(P_2,4);
    r_eff = 2*r_1*r_2/(r_1+r_2);
    h12 = d-(r_1+r_2);
        
    N_overlap = N_overlap+1;
    theta = acosd((d^2+r_1^2-r_2^2)/(2*d*r_1));
    r_c = real(r_1*sind(theta));
    C_List(i,5) = r_c;       
    
    C_c =2*K_s*r_c;
    
    C_s1 = pi*K_s*(Xi*r_eff)^2/r_1;
    C_s2 = pi*K_s*(Xi*r_eff)^2/r_2;
    
    Cond(i) = 1/(1/C_s1+1/C_c+1/C_s2);
    C_List(i,9) = Cond(i);
        
end

Con = sparse(C_List(:,1),C_List(:,2),Cond,N_particles,N_particles);
Con = Con + transpose(Con);

Cona = zeros(N_particles,1);

for i=1:1:N_particles
Cona(i) = sum(Con(i,:));
end

Con = -1*Con;

for i=1:1:N_particles
Con(i,i) = Cona(i);
end

Consave=Con;
% fprintf('Conduction Matrix created\n');

%%
layer_c = 0;

for i = 1:1:N_particles
    if (N(i,direction_cond)-min(N(:,direction_cond)) <= mean_r)%same temperature
        layer_c = layer_c+1;
    end
end
% layer_c 
layer_bot = layer_c;

layer_c = 0;
for i = 1:1:N_particles
    if ((max(N(:,direction_cond))-N(i,direction_cond)) <= mean_r)%same temperature
        layer_c = layer_c+1;
    end
end
% layer_c ;
layer_top = N_particles - layer_c + 1;

Bot = zeros(N_particles,1);
Top = zeros(N_particles,1);
Vol = zeros(N_particles,1);
Inner_part = zeros(N_particles,1);

%%%Voltage initiation
     for i=1:1:layer_bot
        Bot(i,1) = N(i,8);
        Vol(N(i,8))= Top_temp;
     end
    
    for i= layer_top : 1 : N_particles      
         Top(i,1) = N(i,8);
         Vol(N(i,8))=Bot_temp;   
    end
    

    for i= layer_bot+1:1:layer_top-1
        Inner_part(i,1) = N(i,8);
    end
    
Top = nonzeros(Top);
Bot = nonzeros(Bot);
Inner_part = nonzeros(Inner_part);


Y1 = zeros(size(Con,1),1);

    for inner_p = Inner_part
        for tophere=Top
           Y1(inner_p,1) = Y1(inner_p,1)-Con(inner_p,tophere)*Vol(tophere,1);
        end
        for bothere=Bot
           Y1(inner_p,1) = Y1(inner_p,1)-Con(inner_p,bothere)*Vol(bothere,1);
        end        
    end

Bound_here = [Top;Bot]; 

Con_1 = Con;
Con_1(Bound_here,:)=[];
Con_1(:,Bound_here)=[];
Y1(Bound_here,:)=[];

%%
% fprintf('Solving Started\n');
%Solve for unknown Temperatures
X = Con_1\Y1;
% fprintf('Solved\n');

Vol(Inner_part,1)=X;

I1 = zeros(N_particles,1);
I2 = zeros(N_particles,1);

for tophere = Top
   for j=1:1:N_particles
    I1(tophere,1) =  I1(tophere,1) + Con(tophere,N(j,8))*Vol(N(j,8),1);
   end
end

for bothere = Bot
   for j=1:1:N_particles
    I2(bothere,1) =  I2(bothere,1) + Con(bothere,N(j,8))*Vol(N(j,8),1);
   end
end

Avg_current = -(-sum(I1)+sum(I2))/2;

Length_cell =  abs(Wall(direction_cond,2)-Wall(direction_cond,1));
% Area_cell = (max(N(:,1))-min(N(:,1))+2*mean_r)*(max(N(:,2))-min(N(:,2))+2*mean_r);

Length_keff = mean(N(Top,direction_cond))-mean(N(Bot,direction_cond)) ;
Area_cell = abs((Wall(other_dir(1),2)-Wall(other_dir(1),1))*(Wall(other_dir(2),2)-Wall(other_dir(2),1)));

k_eff = (Avg_current*Length_keff)/(Delta_T*Area_cell);

pf = 0;
for i=1:1:N_particles
pf = pf+(4*pi/3)*(N(i,4)^3);
end

pf=pf/(Length_cell*Area_cell);

O_CN = 2*N_overlap/N_particles;
mean_rc = sum(C_List(:,5))/N_overlap;

Ceff = 2*K_s*mean_rc;

k_ana = pf*O_CN*(Ceff)/(2*pi*mean_r);

counter_file =counter_file+1;

File_part(counter_file,1) = counter_file;
File_part(counter_file,2) = direction_cond;
File_part(counter_file,3) = pf;

File_part(counter_file,4) = O_CN;

File_part(counter_file,8) = mean_rc;

File_part(counter_file,10) = mean(C_List(:,4));

File_part(counter_file,11) = sind(mean(C_List(:,4)))^2;


File_part(counter_file,13) = 0;
File_part(counter_file,14) = Length_cell;
File_part(counter_file,15) = Area_cell;
File_part(counter_file,16) = Avg_current;

File_part(counter_file,17) = Ceff*1000;

File_part(counter_file,21) = k_eff;%Numerical Keff
File_part(counter_file,22) = k_ana;%Analytical Keff
File_part(counter_file,23) = (k_eff-k_ana)*100/(k_eff);

fprintf('K_sim \t K_ana \t Error \n');
fprintf('%0.4f\t%0.4f\t%0.4f\n\n',File_part(counter_file,21:23));
File_part(counter_file,24:30)=0;

N_out = N;
N_out(:,5) = Vol;
% N_out(:,1:3) = N_out(:,1:3) + mean_r;
%%%Slopes at each region 25,26,27,28,29
% Average of slopes at 30

%%
%Regional Fitting and conductance
LT=10;
Contact_condi_Lis=zeros(50/LT,10);
C_List2 = C_List;
C_List2(:,1) = C_List(:,2);
C_List2(:,2) = C_List(:,1);
C_List2(:,10) = C_List(:,11);
C_List2(:,11) = C_List(:,10);

C_List3 = [C_List;C_List2];

for il=1:1:50/LT
   
    Fi_1 = N_out(N_out(:,direction_cond) >= ((il-1)*LT*mean_r+Wall(direction_cond,1)),:);
    Fi = Fi_1(Fi_1(:,direction_cond) <= (LT*il*mean_r + Wall(direction_cond,1)),:);
    Fi = Fi(Fi(:,5)>Top_temp-1,:);
    Fit_values = polyfit(Fi(:,direction_cond),Fi(:,5),1);
    File_part(counter_file,24+il) = Avg_current/(Area_cell*Fit_values(1));  

    Fi_2 = C_List3(C_List3(:,10) >= ((il-1)*LT*mean_r+Wall(direction_cond,1)),:);
    Fi_5 = Fi_2(Fi_2(:,10) <= (LT*il*mean_r+Wall(direction_cond,1)),:);

    N_pr = size(Fi,1);
    Contact_condi_Lis(il,1:2) = [il,N_pr];
    Contact_condi_Lis(il,3) = nnz(Fi_5(:,5))/N_pr;
    Contact_condi_Lis(il,6) = sum(Contact_condi_Lis(il,3:5));
    Contact_condi_Lis(il,7) = sum(Fi_5(:,5))/nnz(Fi_5(:,5));
    Contact_condi_Lis(il,10) =  Avg_current/(Area_cell*Fit_values(1));
         

CCr1 = 2*K_s*Contact_condi_Lis(il,6);


Contact_condi_Lis(il,11) = (Contact_condi_Lis(il,3))*( 1/(1/C_s1+1/C_s1+1/CCr1))*1000;

%%
end


File_part(counter_file,30)=mean(File_part(counter_file,25:29));
% File_part = round(File_part*10000)/10000;

% Voltname = sprintf('%s_Kg%0.2f_T_%d.dat',filename,K_s,direction_cond);
% dlmwrite(Voltname,N_out,' ');

% Con_List = sprintf('%s_Kg%0.2f_C_%d.dat',filename,K_s,direction_cond);
% dlmwrite(Con_List,C_List,' ');


end

% Namedat = sprintf('%s_Kg%0.2f_keff.dat',filename,K_s);
% dlmwrite(Namedat,File_part,' ');


toc;
