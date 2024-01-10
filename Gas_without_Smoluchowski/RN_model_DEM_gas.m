clear all;
tic;

% fileID = fopen('S50P.dat','r');
%  [ScFA,ScFCount] = fscanf(fileID,'%s',7);
% fclose(fileID);
% double n;
% ScF = str2double(ScFA(30:36));
Ann_File = zeros(1,10);
%%
Tc = 25;
Tk = Tc+273.14;
%Li4OSi4 
% K_s = (((7.317E-12)*(Tc.^(4.0)))+((-1.302E-8)*(Tc.^(3.0)))+(8.712E-6*(Tc.^(2.0)))+(-0.002876*Tc)+2.620);

%Li2TiO3
poro = 0.1;
kbeta = 1.06-(2.88E-4)*Tk;
K_s = ((1-poro)/(1+kbeta*poro))*(4.77 - (5.11E-3)*Tk + (3.12E-6)*Tk*Tk);

%Helium
K_f = 3.366*(10^(-3))*(Tk.^(0.668));

%Air
% K_f =(-1*10^(-11))*(Tc.^3)  - (4*(10^(-8))*(Tc.^2)) + (8*(10^(-5))*Tc) + 0.0241;

%%
Start = 80;
Jump = 5;
End = 80;
Number_files = floor((End-Start)/Jump +1);

direction_cond = 3;

Bthres=100;
L_thres = 1;
Xi = 0.71;
cutoff_factor = 0.5;
 
requiredRk = 1E-3; %Scaling doesnt change anything--for gas case 
%without smoluchowski effect and DEM is independent of radius.
%keff independent of radius of pebbles.

Top_temp = 10;
Bot_temp = 20;

D_Strain = 0.00015;
alpha = K_s/K_f;

Delta_T = abs(Top_temp-Bot_temp);
     
File_part = zeros(Number_files,20);
counter_file =0;
       
for file_number = Start:Jump:End

    fprintf('file number = %d\n\n',file_number)
    


 %%
%Reading the file
   if file_number<10
    file_name = 'S0%dP.dat';
    else 
    file_name = 'S%dP.dat';
    end

file_name_str = sprintf(file_name,file_number);
N = dlmread(file_name_str,' ',2,0);
datatxt = dlmread(file_name_str,' ',[1 1 2 11]);
% fileID1 = fopen(file_name_str,'r');
% [datatxt,ScFCount] = fscanf(fileID1,'%s',50);
% fclose(fileID1);
Straink = datatxt(1,3)*-1;
Stressk = datatxt(1,6)/-1E6;
Rmaxk = max(N(:,4));

Scfk = requiredRk/Rmaxk;
ScF = Scfk;


% ScF=0.016067;
% Scfk=ScF;
%Getting the scaling factor


%Adjusting with Scaling Factor
N(:,1:4) = Scfk*N(:,1:4);

%Sorting through the columns based on direction
N = sortrows(N,direction_cond);

N_particles = size(N,1);
mean_r  = mean(N(:,4));

N(:,8) = 1:1:N_particles;

% fprintf('Particle DATA is read\n');
%%
%Finding the periodic boundary particles

Other_dir = [1;2;3];
Other_dir(direction_cond) = [];

F1 = Other_dir(1);R1 = Other_dir(2);

PF = N(N(:,F1) < min(N(:,F1)) + 1.1*mean_r,:);
PB = N(N(:,F1) > max(N(:,F1)) - 1.1*mean_r,:);
PR = N(N(:,R1) < min(N(:,R1)) + 1.1*mean_r,:);
PL = N(N(:,R1) > max(N(:,R1)) - 1.1*mean_r,:);

PF_per = PF;
PB_per = PB;
PR_per = PR;
PL_per = PL;

PF_per(:,F1) = PF_per(:,F1) + ScF;
PB_per(:,F1) = PB_per(:,F1) - ScF;
PR_per(:,R1) = PR_per(:,R1) + ScF;
PL_per(:,R1) = PL_per(:,R1) - ScF;


%%
%Creating the contact list
M1 = contact_list(N,cutoff_factor,direction_cond);
M2 = contact_list_B(PF,PB_per,cutoff_factor);
M3 = contact_list_B(PR,PL_per,cutoff_factor);

C_List = [M1;M2;M3];
% C_List = M1;

N_contacts = size(C_List,1);

% fprintf('Contact list is made\n');

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
    
    if h12<0
        
        N_overlap = N_overlap+1;
        theta = acosd((d^2+r_1^2-r_2^2)/(2*d*r_1));
        r_c = real(r_1*sind(theta));
        C_List(i,5) = r_c;
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
    if (N(i,direction_cond) <= mean_r)%same temperature
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

Avg_current = (-sum(I1)+sum(I2))/2;

Length_cell =  ScF*(1-Straink);
Area_cell = ScF^2;
Length_keff = mean(N(Top,direction_cond))-mean(N(Bot,direction_cond));
k_eff = -(Avg_current*Length_keff)/(Delta_T*Area_cell);

pf = 0;
for i=1:1:N_particles
pf = pf+(4*pi/3)*(N(i,4)^3);
end

pf=pf/(Length_cell*Area_cell);

if(file_number==Start)
    pforgk = pf;
end

O_CN = 2*N_overlap/N_particles;
G_CN = 2*N_gap/N_particles;
T_CN = 2*N_touch/N_particles;

CN_t = O_CN+G_CN+T_CN;

mean_rc = sum(C_List(:,5))/N_overlap;
Eps = sum(C_List(:,8))/N_gap;

eff_h = (Xi^2)*(r_eff)/(exp(Eps)-1);

     beta = mean_rc*alpha/mean_r;        
         
        if beta < 1
          Kc = 0.22*(beta^2); DKg = -0.05*(beta^2);
        else if beta < Bthres
               Kc =  (beta-1)*(2*Bthres/pi-0.22)/(Bthres-1)+0.22 ;DKg =(beta-1)*(-2*log(Bthres)+0.05)/(Bthres-1)-0.05;
            else 
             Kc = 2*(beta)/pi; DKg = -2*log(beta);
            end
        end
CC1 = pi*K_f*mean_r*(Kc+DKg+log(alpha^2));
CC3 = pi*K_f*mean_r*(log(alpha^2)+log(1+Xi^2*alpha^2))/2;
CC2 = pi*K_f*mean_r*Eps; 
C_s1 =  pi*K_s*(Xi*mean_r)^2/mean_r;
Ceff = (O_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC1)) +  ((G_CN+T_CN)/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC2)) ; 

k_ana = pf*CN_t*(Ceff)/(2*pi*mean_r);
counter_file = counter_file+1;

Ann_File(end+1,1) = pforgk;
Ann_File(end,2) = Stressk;
Ann_File(end,3) = Straink;
Ann_File(end,4) = requiredRk;
Ann_File(end,5) = K_s;
Ann_File(end,6) = K_f;
Ann_File(end,7) = k_eff;
Ann_File(end,8) = k_ana;
Ann_File(end,9) = pf;
Ann_File(end,10) = file_number;
File_part(counter_file,1) = file_number;
File_part(counter_file,2) = file_number*D_Strain*100;
File_part(counter_file,3) = pf;

File_part(counter_file,4) = CN_t;
File_part(counter_file,5) = O_CN/CN_t;
File_part(counter_file,6) = G_CN/CN_t;
File_part(counter_file,7) = T_CN/CN_t;

File_part(counter_file,8) = mean_rc/mean_r;
File_part(counter_file,9) = Eps;
File_part(counter_file,10) = mean(C_List(:,4));

File_part(counter_file,11) = sind(mean(C_List(:,4)))^2;
File_part(counter_file,12) = eff_h/r_eff;
File_part(counter_file,13) = k_eff;
File_part(counter_file,14) = k_ana;
File_part(counter_file,15) = (k_eff-k_ana)*100/(k_eff);

File_part(counter_file,16) = Ceff*1000;
File_part(counter_file,17) = 0;
File_part(counter_file,18) = (O_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC1))*1000;
File_part(counter_file,19) = (G_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC2))*1000;
File_part(counter_file,20) = (T_CN/CN_t)*( 1/(1/C_s1+1/C_s1+1/CC3))*1000;
Tco = File_part(counter_file,4).*(File_part(counter_file,18)+File_part(counter_file,19)+File_part(counter_file,20));
fprintf('K_sim \t K_ana \t Error\n');
fprintf('%0.4f\t%0.4f\t%0.4f\n\n',File_part(counter_file,13:15));

  if file_number<10
    formatSpec3 = 'H=%.3f_T0%dCR.dat';
  
    else 
    formatSpec3 = 'H=%.3f_T%dCR.dat';
    
  end
    

str1 = sprintf(formatSpec3,cutoff_factor,file_number);
% dlmwrite(str1,C_List,' ');

end

inipf=(5000*4*pi*(mean_r)^3)/(3*ScF^3);
Namedat = sprintf('Kg%0.2f_Kf%0.3f_pf%.3f_Keff.dat',K_s,K_f,inipf);
% dlmwrite(Namedat,File_part,' ');


% save('Ann_1.mat','Ann_File');

toc;