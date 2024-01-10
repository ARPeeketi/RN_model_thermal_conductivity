clear all;
tic;

%%
Pressure = 1E5;
Tc = 25;
Tk = Tc+273.14;

       p1 =  -1.596e-09 ;
       p2 =   3.042e-06 ;
       p3 =   -0.001923  ;
       p4 =       2.591 ;
% 
% %Li4OSi4 
K_s = p4 + p3*Tc + p2*Tc^2 + p1*Tc^3;
mg = 119.85;
K_s=150;

 Xi = 0.36;
% % Li2TiO3
% 
% poro = 0.1;
% kbeta = 1.06-(2.88E-4)*Tk;
% K_s = ((1-poro)/(1+kbeta*poro))*(4.77 - (5.11E-3)*Tk + (3.12E-6)*Tk*Tk);
% mg=109.76;


%Helium
K_f0 = 3.366*(10^(-3))*(Tk.^(0.668));
mf=4;
dm = 2.15E-10;
visc = 1.865*(1e-5)*(Tk/273.16)^(0.7);


mr=mg/mf;
ac=2.4*(mr)/((1+mr)^2);
bsbeta=(19/12)*(2-ac)/ac;

kBo = 1.38E-23;

 %%
filename = 'FFSN1_05'; %%%%Change file name when required...
fileread = sprintf('%s.dat',filename);
N = dlmread(fileread);
ScF=1;

N(:,1:4) = ScF*N(:,1:4);

% [a0w,b0w] = max(N(:,3));
% [a0m,b0m] = min(N(:,3));
a0m = 5.0152E-4;
a0w = 0.0267;
Wall=[-0.0125,0.0125;-0.0125,0.0125;a0m-N(1,4),a0w+N(1,4)];

fprintf('Particle DATA is read\n\n');


cutoff_factor = 0.5;
Bthres=100;
L_thres = 1;

Top_temp = 10;
Bot_temp = 20;

Delta_T = abs(Top_temp-Bot_temp);
counter_file = 0;

dir_other = [2,3;1,3;1,2];
File_part = zeros(1,20);

%%
%Change the below for loop for specific direction
for direction_cond = 3:1:3

other_dir = dir_other(direction_cond,:);   


%Sorting through the columns based on direction
N = sortrows(N,direction_cond);
disp(unique(N(:,4)));
N_particles = size(N,1);
mean_r  = mean(N(:,4));

N(:,8) = 1:1:N_particles;

%%
%Creating the contact list
C_List = contact_list(N,cutoff_factor,direction_cond);

N_contacts = size(C_List,1);

fprintf('Contact list is generated for direction %d\n',direction_cond);

Cond = zeros(N_contacts,1);

N_overlap=0;
N_gap=0;
N_touch = 0;

%Conduction Matrix

for i=1:1:N_contacts
    
    P_1 = C_List(i,1); 
    P_2 = C_List(i,2);   
    d = C_List(i,3);
    
    r_1 = N(P_1,4);
    r_2 = N(P_2,4);
    r_eff = 2*r_1*r_2/(r_1+r_2);
    h12 = d-(r_1+r_2);
    rmin = min(r_1,r_2);
    rmax = min(r_1,r_2);
    
   if (Xi*r_eff < rmin)
    Gamma_1 = asind((Xi*r_eff)/r_1); 
    Gamma_2 = asind((Xi*r_eff)/r_2);
   else
       Gamma_1 = asind(rmin/r_1);
       Gamma_2 = asind(rmin/r_2);
   end
   
    if h12<0

        N_overlap = N_overlap+1;
        
        r_c = sqrt(abs(h12)*r_eff/2);
        
        theta1 = asind(r_c/r_1);
        theta2 = asind(r_c/r_2);
        C_List(i,5) = r_c;
                
        Lo1 = r_1*( 1- 0.5*( cos(theta1) + cos(Gamma_1) ) );
        Lo2 = r_2*( 1- 0.5*( cos(theta2) + cos(Gamma_2) ) );
%         Lo1 = r_1*(1- (4/3)*( ((cos(theta1))^3 - (cos(Gamma_1))^3)/(cos(2*theta1) - cos(2*Gamma_1)) ));
%         Lo2 = r_1*(1- (4/3)*( ((cos(theta2))^3 - (cos(Gamma_2))^3)/(cos(2*theta2) - cos(2*Gamma_2)) ));
        
        Lo=Lo1+Lo2 + h12;
        
        
        Hff = (visc/(Pressure*Lo))*sqrt(pi*8314*Tk/(2*mf));
        K_f = K_f0/(1+(2*bsbeta*Hff));
        C_List(i,12) = Lo;
        alpha = K_s/K_f;
        beta = alpha*r_c/r_eff;
        
        if beta < 1
          Kc = 0.22*(beta^2); DKg = -0.05*(beta^2);
        else if beta < Bthres
               Kc =  (beta-1)*(2*Bthres/pi-0.22)/(Bthres-1)+0.22 ;DKg =(beta-1)*(-2*log(Bthres)+0.05)/(Bthres-1)-0.05;
            else 
             Kc = 2*(beta)/pi; DKg = -2*log(beta);
            end
        end
        
        C_c = pi*K_f*r_eff*(Kc+DKg+log(alpha^2));
    else
        Lg1 = r_1*( 1- 0.5*( 1 + cos(Gamma_1) ) );
        Lg2 = r_2*( 1- 0.5*( 1 + cos(Gamma_2) ) );
%         Lg1 = r_1*(1 - (4/3)*(((1-(cos(Gamma_1))^3))/(1-cos(2*Gamma_1))));
%         Lg2 = r_2*(1 - (4/3)*(((1-(cos(Gamma_2))^3))/(1-cos(2*Gamma_2))));
        
        Lg = Lg1+Lg2+h12;
        C_List(i,13) = Lg;
        Hff = (visc/(Pressure*Lg))*sqrt(pi*8314*Tk/(2*mf));
        K_f = K_f0/(1+(2*bsbeta*Hff));
        alpha = K_s/K_f;
        lamda = (alpha^2)*h12/r_eff;
        
        if lamda < L_thres
            C_List(i,6) = h12;
            N_touch = N_touch+1;
%             C_c = pi*K_f*r_eff*log(alpha^2);
            C_c = pi*K_f*r_eff*(lamda*(log(1+Xi^2*alpha^2)-log(alpha^2)) + log(alpha^2));
        else
            N_gap = N_gap+1;
            C_List(i,7) = h12;
            C_List(i,8) = log(1 + (Xi^2)*r_eff/h12);
            C_c = pi*K_f*r_eff*log(1 + (Xi^2)*r_eff/h12);
        end
    end
    
    C_s1 = pi*K_s*(Xi*r_eff)^2/r_1;
    C_s2 = pi*K_s*(Xi*r_eff)^2/r_2;
    
    Cond(i) = 1/(1/C_s1+1/C_c+1/C_s2);
    C_List(i,9) = Cond(i);
        
%     if (mod(i,5000)==0)
%     disp(i);
%     end
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
G_CN = 2*N_gap/N_particles;
T_CN = 2*N_touch/N_particles;

CN_t = O_CN+G_CN+T_CN;
mean_rc = sum(C_List(:,5))/N_overlap;
[minrc,minrcid] = min(abs(C_List(:,5)-mean_rc));
Loo = C_List(minrcid,12);

Hff = (visc/(Pressure*Loo))*sqrt(pi*8314*Tk/(2*mf));
K_fo = K_f0/(1+(2*bsbeta*Hff));

Eps = sum(C_List(:,8))/N_gap;
eff_h = (Xi^2)*(r_eff)/(exp(Eps)-1);
[minhe,minheid] = min(abs(C_List(:,7)-eff_h));
Lgg = C_List(minheid,13);

Hff = (visc/(Pressure*Lgg))*sqrt(pi*8314*Tk/(2*mf));
K_fg = K_f0/(1+(2*bsbeta*Hff));

    alpha = K_s/K_fo;
    beta = mean_rc*alpha/mean_r;        
         
    if beta < 1
      Kc = 0.22*(beta^2); DKg = -0.05*(beta^2);
    else if beta < Bthres
           Kc =  (beta-1)*(2*Bthres/pi-0.22)/(Bthres-1)+0.22 ;DKg =(beta-1)*(-2*log(Bthres)+0.05)/(Bthres-1)-0.05;
        else 
         Kc = 2*(beta)/pi; DKg = -2*log(beta);
        end
    end
        
CC1 = pi*K_fo*mean_r*(Kc+DKg+log(alpha^2));
alpha = K_s/( (K_fo+K_fg)/2 );
CC3 = pi*((K_fo+K_fg)/2)*mean_r*(log(alpha^2)+log(1+Xi^2*alpha^2))/2;
CC2 = pi*K_fg*mean_r*Eps; 

Ceff = (O_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC1)) +  ((G_CN)/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC2)) + (T_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC3)) ; 

k_ana = pf*CN_t*(Ceff)/(2*pi*mean_r);
counter_file =counter_file+1;

File_part(counter_file,1) = counter_file;
File_part(counter_file,2) = direction_cond;
File_part(counter_file,3) = pf;

File_part(counter_file,4) = CN_t;
File_part(counter_file,5) = O_CN/CN_t;
File_part(counter_file,6) = G_CN/CN_t;
File_part(counter_file,7) = T_CN/CN_t;

File_part(counter_file,8) = mean_rc/mean_r;
File_part(counter_file,9) = Loo/mean_r;
File_part(counter_file,10) = Eps;

File_part(counter_file,11) = eff_h/mean_r;
File_part(counter_file,12) = Lgg/mean_r;

File_part(counter_file,13) = 0;
File_part(counter_file,14) = Length_cell;
File_part(counter_file,15) = Area_cell;
File_part(counter_file,16) = Avg_current;

File_part(counter_file,17) = Ceff*1000;
File_part(counter_file,18) = (O_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC1))*1000;
File_part(counter_file,19) = (G_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC2))*1000;
File_part(counter_file,20) = (T_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC3))*1000;

File_part(counter_file,21) = k_eff;%Numerical Keff
File_part(counter_file,22) = k_ana;%Analytical Keff
File_part(counter_file,23) = (k_eff-k_ana)*100/(k_eff);

fprintf('K_sim \t K_ana \t Error \n');
fprintf('%0.4f\t%0.4f\t%0.4f\n\n',File_part(counter_file,21:23));
File_part(counter_file,24:30)=0;

N_out = N;
N_out(:,5) = Vol;
plot(N_out(:,3),N_out(:,5),'r.');
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
    Fit_values = polyfit(Fi(:,direction_cond),Fi(:,5),1);
    File_part(counter_file,24+il) = Avg_current/(Area_cell*Fit_values(1));  

%     Fi_2 = C_List3(C_List3(:,10) >= ((il-1)*LT*mean_r+Wall(direction_cond,1)),:);
%     Fi_5 = Fi_2(Fi_2(:,10) <= (LT*il*mean_r+Wall(direction_cond,1)),:);
% 
%     N_pr = size(Fi,1);
%     Contact_condi_Lis(il,1:2) = [il,N_pr];
%     Contact_condi_Lis(il,3:5) = [nnz(Fi_5(:,5))/N_pr,nnz(Fi_5(:,6))/N_pr,nnz(Fi_5(:,7))/N_pr];
%     Contact_condi_Lis(il,6) = sum(Contact_condi_Lis(il,3:5));
%     Contact_condi_Lis(il,7) = sum(Fi_5(:,5))/nnz(Fi_5(:,5));
%     Contact_condi_Lis(il,8) = sum(Fi_5(:,8))/nnz(Fi_5(:,8)); 
%     Contact_condi_Lis(il,9) = (Xi^2)*(mean_r)/(exp(Contact_condi_Lis(il,8))-1);
%     Contact_condi_Lis(il,10) =  Avg_current/(Area_cell*Fit_values(1));
%     
%     Rmean_rc = Contact_condi_Lis(il,7);
%     [Rminrc,Rminrcid] = min(abs(Fi_5(:,5)-Rmean_rc));
%     RLoo = Fi_5(Rminrcid,12);
% 
%     Hff = (visc/(Pressure*RLoo))*sqrt(pi*8314*Tk/(2*mf));
%     K_for = K_f0/(1+(2*bsbeta*Hff));
% 
%     REps = Contact_condi_Lis(il,8);
%     Reff_h = Contact_condi_Lis(il,9);
%     [Rminhe,Rminheid] = min(abs(Fi_5(:,7)-Reff_h));
%     RLgg = Fi_5(Rminheid,13);
%     Hff = (visc/(Pressure*RLgg))*sqrt(pi*8314*Tk/(2*mf));
%     K_fgr = K_f0/(1+(2*bsbeta*Hff));
% 
%     alpha = K_s/K_for;
%          beta = Contact_condi_Lis(il,7)*alpha/mean_r;        
%          
%         if beta < 1
%           Kc = 0.22*(beta^2); DKg = -0.05*(beta^2);
%         else if beta < Bthres
%                Kc =  (beta-1)*(2*Bthres/pi-0.22)/(Bthres-1)+0.22 ;DKg =(beta-1)*(-2*log(Bthres)+0.05)/(Bthres-1)-0.05;
%             else 
%              Kc = 2*(beta)/pi; DKg = -2*log(beta);
%             end
%         end
% CCr1 = pi*K_for*mean_r*(Kc+DKg+log(alpha^2));
% alpha = K_s/((K_for+K_fgr)/2);
% CCr2 = pi*((K_for+K_fgr)/2)*mean_r*(log(alpha^2)+log(1+Xi^2*alpha^2))/2;
% CCr3 = pi*K_fgr*mean_r*Contact_condi_Lis(il,8); 
% 
% Contact_condi_Lis(il,11) = (Contact_condi_Lis(il,3))*( 1/(1/C_s1+1/C_s1+1/CCr1))*1000;
% Contact_condi_Lis(il,12) = (Contact_condi_Lis(il,4))*( 1/(1/C_s1+1/C_s1+1/CCr2))*1000;
% Contact_condi_Lis(il,13) = (Contact_condi_Lis(il,5))*( 1/(1/C_s1+1/C_s1+1/CCr3))*1000;
%%
end


File_part(counter_file,30)=mean(File_part(counter_file,25:29));

Voltname = sprintf('%s_Kg%0.2f_Kf%0.2f_T_%d.dat',filename,K_s,K_f,direction_cond);
% dlmwrite(Voltname,N_out,' ');

Con_List = sprintf('%s_Kg%0.2f_Kf%0.2f_C_%d.dat',filename,K_s,K_f,direction_cond);
% dlmwrite(Con_List,C_List,' ');
end

File_part = round(File_part*10000)/10000;
Namedat = sprintf('%s_Kg%0.2f_Kf%0.2f_1.dat',filename,K_s,K_f);
% dlmwrite(Namedat,File_part,' ');

toc;