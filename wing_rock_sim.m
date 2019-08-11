clc;
Kp=1.5;
Kd=1.3;
ti=0;
tf=60;
dt=0.005;
Am=[[0 1];[-Kp -Kd]];
% First section of the code simply demonstates PD controller response
% to step inputs. This is the ref model that we wish to follow 

phi_p_prev=[0;0]; % Initialised vector
phi_p_ref_store=[]; % matrix for storing phi and phidot (p) values as vector columns
phi_des_store=[];

for i=ti:dt:tf
    if i>10 && i<20  % The desired signal we wish to track
        phi_des=1;
        phid_des=0;
        phidd_des=0;
    elseif i>30 && i<40
        phi_des=-1;
        phid_des=0;
        phidd_des=0;
    else
        phi_des=0;
        phid_des=0;
        phidd_des=0;
    end
            
    
    phi_des_store=[phi_des_store,phi_des];
    
    phi_p_ref_next=phi_p_prev+dt*(Am*phi_p_prev+[0;1]*...
    (phidd_des+Kd*phid_des+Kp*phi_des));
    
    phi_p_ref_store=[phi_p_ref_store,phi_p_prev];
    phi_p_prev=phi_p_ref_next;
end

% for plotting reference signal, PD ref model response

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_vec=ti:dt:tf;
figure(1)
plot(time_vec,phi_des_store,'r--');
hold on;
grid on;
plot(time_vec,phi_p_ref_store(1,:));
xlabel('Time')
ylabel('ref signal and PD response')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Actual model with adaptive control (using above control law as linear model)

phi_p_prev_2=[0;0]; % Initialised vector
phi_p_adaptive_store=[]; 
phi_des_store_2=[];

BW=2;
centers=10;
gain=200;
arbit_centres=repmat(linspace(-0.05,0.05,centers),[2,1]); 
% repmat(linspace(3,5,10),[2,1]);  

% center location change causes issue. So if domain of operation uncertain,
% issue arises. Shortcoming of classical MRAC can be seen by replacing
% value of arbit_centres with the above commented repmat command

% Also making vad=0, will show how badly PD response is when delta is
% present

basis=zeros(centers,1);
counter=0;
outputweight=zeros(centers,1);
B=[0;1];
Q=eye(2);
P=lyap(Am',Q); % Am'P+P*Am+Q = 0
delta_store=[]; % for plotting delta uncertainity
adaptive_store=[]; % for plotting adaptive element

for i=ti:dt:tf
    counter=counter+1;
     if i>10 && i<20
        phi_des_2=1;
        phid_des_2=0;
        phidd_des_2=0;
    elseif i>30 && i<40
        phi_des_2=-1;
        phid_des_2=0;
        phidd_des_2=0;
    else
        phi_des_2=0;
        phid_des_2=0;
        phidd_des_2=0;
     end
    
    phi_des_store_2=[phi_des_store_2,phi_des_2];
    %%%%%%%%%%%% adaptive element %%%%%%%%%%%%%
    for m=1:centers
        basis(m)=exp((-norm(arbit_centres(:,m)-phi_p_prev_2)^2)/(2*BW));
    end
    error=phi_p_ref_store(:,counter)-phi_p_prev_2;
    
    outputweight=outputweight+dt*(-gain)*basis*error'*P*B;
    
    vad=outputweight'*basis;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Uncertainity
    delta=0.8+0.2314*phi_p_prev_2(1)+0.6918*phi_p_prev_2(2)-0.6245*abs(phi_p_prev_2(1))*phi_p_prev_2(2)...
        +0.0095*abs(phi_p_prev_2(2))*phi_p_prev_2(2)+0.0214*(phi_p_prev_2(1))^3;
    
    % Adaptive element countering uncertainity
    phi_p_next_2=phi_p_prev_2+dt*([[0 1];[-Kp -Kd]]*phi_p_prev_2+[0;1]*...
    (phidd_des_2+Kd*phid_des_2+Kp*phi_des_2)+[0;delta]+[0;-vad]);
    
    
    %%%%%%%%%% For plotting later on %%%%%%%%%%
    phi_p_adaptive_store=[phi_p_adaptive_store,phi_p_prev_2];
    delta_store=[delta_store,delta];
    adaptive_store=[adaptive_store,vad];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    phi_p_prev_2=phi_p_next_2;
end

figure(2)
time_vec=ti:dt:tf;
plot(time_vec,phi_des_store_2,'r--','LineWidth',2);
hold on
grid on;
plot(time_vec,phi_p_adaptive_store(1,:));
xlabel('time');
legend('Phi desired','Phi after countering delta with adaptive element')


figure(3)
plot(time_vec,delta_store,'r--','LineWidth',2);
hold on
grid on;
plot(time_vec,adaptive_store);
xlabel('time');
legend('delta','adaptive control element response')
