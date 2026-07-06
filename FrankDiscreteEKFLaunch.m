%% 
clearvars;
close('all');
% % Input Data
InputData = readtable("Lot24CRCGPS - Sheet1.csv");
Lat = table2array(InputData(:,1));
Lon = table2array(InputData(:,2));
XdSint = table2array(InputData(:,3));
YdSint = table2array(InputData(:,4));
XcovInt = table2array(InputData(:,5));
YcovInt = table2array(InputData(:,6));
utmzone = 16;
[xint,yint,utmzone] = deg2utm(Lat(1:length(Lat),:),Lon(1:length(Lon),:));
x = -(yint-yint(1,1));
y = (xint-xint(1,1));
inc = 0;

%ID = readtimetable('LatRun2Velocity2')
InputData = readtable("Lot24CRCINS - Sheet1.csv");
legend("color1","LineColor")

%Xquat = smoothdata((InputData(:,5)), "movmean", 200);
Xquat = smoothdata(table2array(InputData(:,1)), "movmean", .5);
Yquat = smoothdata(table2array(InputData(:,2)), "movmean", .5);
Zquat = smoothdata(table2array(InputData(:,3)), "movmean", .5);
Wquat = smoothdata(table2array(InputData(:,4)), "movmean", .5);
PsidSint = smoothdata(table2array(InputData(:,5)), "movmean", 100);
LatAccelInt = smoothdata(table2array(InputData(:,7)), "movmean", 100);
LongAccelInt = smoothdata(table2array(InputData(:,6)), "movmean", 100);



InputData = readtable("GBCRC1550ENC - Sheet1.csv");
SAint = smoothdata(table2array(InputData(:,6)), "movmean", .5);
TicksYRint = smoothdata(table2array(InputData(:,9)), "movmean", 50);
%Psi =  smoothdata(table2array(InputData(:,1)), "movmean", 200);
% figure
% plot(Xquat)
% hold on
% plot(Yquat)
% plot(Zquat)
% plot(Wquat)
% legend
% hold off
%plot(Psi)

size1 = size(Xquat);
quatArray = zeros(size1(1),3)

% for n = 1:size(Xquat)
%     quatArray(n,1) = Xquat(n,1);
%     quatArray(n,2) = Yquat(n,1);
%     quatArray(n,3) = Zquat(n,1);
%     quatArray(n,4) = Wquat(n,1);
% end
sizexquat = size(Xquat(:,1))
for n = 1:sizexquat(1)
    quat = [Wquat(n,1),Xquat(n,1),Yquat(n,1),Zquat(n,1)];
    quatArray(n,:) = quat2eul(quat,"ZYX");
end
sizequatarray = size(quatArray(:,1))
for n = 2:sizequatarray(1)
    if quatArray(n-1,1)-quatArray(n,1) > 6 %% +->-
        inc = inc + 2*pi;
        Psi(n,:) = quatArray(n,1)+inc;
    elseif quatArray(n-1,1)-quatArray(n,1) < -6 %% -->+
        inc = inc-2*pi;
        Psi(n,:) = quatArray(n,1)+inc;
    else 6 >= quatArray(n-1,1)-quatArray(n,1) >= 6 
        Psi(n,:) = quatArray(n,1)+inc;
    end
end
%plot(Psi)
% figure
% plot(quatArray(:,1))
% hold on
% plot(quatArray(:,2))
% plot(quatArray(:,3))
% plot(PsidSint)
% legend
% hold off 
%plot(quatArray(:,1))

PsiS = timetable(Psi, 'TimeStep', seconds(.01));
plot(Psi)
LatAccel = timetable(LatAccelInt(:,1), 'TimeStep', seconds(.01));
LongAccel = timetable(LongAccelInt(:,1), 'TimeStep', seconds(.01));
PsidS = timetable(PsidSint(:,1), 'TimeStep', seconds(.01));
Pitch = timetable(quatArray(:,2), 'TimeStep', seconds(.01));
Roll = timetable(quatArray(:,3), 'TimeStep', seconds(.01));
SteeringAngle = timetable(SAint(:,1), 'TimeStep', seconds(.01));
TicksYR = timetable(TicksYRint(:,1), 'TimeStep', seconds(.01));
XdS = timetable(XdSint(:,1), 'TimeStep', seconds(.107));
YdS = timetable(YdSint(:,1), 'TimeStep', seconds(.107));
Xcov = timetable(XcovInt(:,1), 'TimeStep', seconds(.107));
Ycov = timetable(YcovInt(:,1), 'TimeStep', seconds(.107));
X = timetable(x, 'TimeStep', seconds(.05));
Y = timetable(y, 'TimeStep', seconds(.05));
figure
plot(X.Variables,Y.Variables)
title('Steady State Cornering Maneuver Path')
xlabel('X (m)');
ylabel('Y (m)');




figure
%plot(PsiS.Time, PsiS.Variables*.1)
%hold on
% plot(XdS.Time, XdS.Variables);
% plot(YdS.Time, YdS.Variables);
plot(Pitch.Time, Pitch.Variables)
 hold on
plot(Roll.Time, Roll.Variables)
% plot(PsidS.Time, PsidS.Variables)
%plot(SteeringAngle.Time, SteeringAngle.Variables)
legend('Yaw','Pitch','roll','Angular Velocity', 'Steering Angle')
hold off 

out = sim("FrankSensingDiscreteEKF.slx")

FfrInt = get(out.logsout,45);
FflInt = get(out.logsout,44);
FrrInt = get(out.logsout,47);
FrlInt = get(out.logsout,46);

TotalForceInt = get(out.logsout,48);
CurveFitForceInt = get(out.logsout,18);
TotalVelInt = get(out.logsout,39);


Alpha2Int = get(out.logsout,32);
Alpha3Int = get(out.logsout,34);
Alpha4Int = get(out.logsout,36);
Alpha5Int = get(out.logsout,37);
AlphaCInt = get(out.logsout,38);

Force1int = get(out.logsout,9);
Force2int = get(out.logsout,10);
Force3int = get(out.logsout,11);
Force4int = get(out.logsout,12);
Force5int = get(out.logsout,13);

Slipint = get(out.logsout,19);

Ffr = FfrInt.Values.Data;
Ffl = FflInt.Values.Data;
Frr = FrrInt.Values.Data;
Frl = FrrInt.Values.Data;

TotalForce = TotalForceInt.Values.Data;
TotalVel = TotalVelInt.Values.Data;
CurveFitForce = CurveFitForceInt.Values.Data;


Alpha2 = Alpha2Int.Values.Data;
Alpha3 = Alpha3Int.Values.Data;
Alpha4 = Alpha4Int.Values.Data;
Alpha5 = Alpha5Int.Values.Data;
AlphaC = AlphaCInt.Values.Data;

Force1 = Force1int.Values.Data;
Force2 = Force2int.Values.Data;
Force3 = Force3int.Values.Data;
Force4 = Force4int.Values.Data;
Force5 = Force5int.Values.Data;

Slip = Slipint.Values.Data;

AlphaTime = Alpha5Int.Values.Time;

% figure
% plot(AlphaTime,TotalVel)
% title('Slip Sweep Velocity Profile')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% xlim([75 150])
% % ylim([0 40])
% grid on
% plot([0 0], ylim, '--')                 % Dashed Vertical Line at x=0
% plot(xlim, [0 0], '--')                 % Dashed Horizontal Line at y=0
% 
% figure
% % plot(FflInt.Values.Time,FflInt.Values.Data);
% plot(AlphaTime(73001:150001),Alpha2(73001:150001));
% hold on
% plot(FfrInt.Values.Time,FfrInt.Values.Data);
% plot(FrrInt.Values.Time,FrrInt.Values.Data);
% plot(FrlInt.Values.Time,FrlInt.Values.Data);
figure
plot(AlphaTime(73001:150001),AlphaC(73001:150001),'LineWidth',1,'Color',"b",'LineStyle','--');
hold on
plot(AlphaTime(73001:150001),Alpha2(73001:150001));
plot(AlphaTime(73001:150001),Alpha3(73001:150001));
plot(AlphaTime(73001:150001),Alpha4(73001:150001));
plot(AlphaTime(73001:150001),Alpha5(73001:150001));
plot(xlim, [0 0], '--', 'HandleVisibility','off','Color',[0 0 0])                 % Dashed Horizontal Line at y=0
grid on
xlim([75 150])
ylim([-.05 .3])
legend('Center of Mass','Wheel 2','Wheel 3','Wheel 4','Wheel 5',' ')
title('Lateral Wheel Slip over Slip Sweep Maneuver')
xlabel('Time (s)')
ylabel('Wheel Slip (rad)')
hold off

figure
plot(Slip,Force1,'LineWidth',2)
hold on
plot(Slip,Force2,'LineWidth',2)
plot(Slip,Force3,'LineWidth',2)
plot(Slip,Force4,'LineWidth',2)
plot(Slip,Force5,'LineWidth',2)
xlim([0 1])
ylim([-10 70])
title('Pure Lateral Slip-Force Curves')
xlabel('Slip (rad)')
ylabel('Tire Lateral Force (N)')
legend('Fz = 50N','Fz = 45N','Fz = 40N','Fz = 35N','Fz = 30N')
hold off

figure

% plot(AlphaC(83001:150001), TotalForce(83001:150001)*.25);
% title('Steady State Cornering Slip and Tire Force')
% xlabel('Wheel Slip (rad)')
% ylabel('Tire Force (N)')
% xlim([-.05 .2])
% ylim([0 40])
% grid on
% hold on
% plot([0 0], ylim, '--')                 % Dashed Vertical Line at x=0
% plot(xlim, [0 0], '--')                 % Dashed Horizontal Line at y=0
% hold off



xdata = ...
 [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = ...
 [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];

fun = @(x,xdata)x(1)*exp(x(2)*xdata);
x0 = [100,-1];
x = lsqcurvefit(fun,x0,xdata,ydata)



DataX = [Alpha2,Alpha3,Alpha4,Alpha5,Ffr,Ffl,Frr,Frl,AlphaC];

% C = 1;
% D = 2*Ffr.^2+3*Ffr;    
% BCD = 4*sin(5*atan(6*Ffr));
% B = (4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1);
% E = 7*Ffr.^2+8*Ffr+9;
mI = 73001;
mF = 150001;

%Ffr = (1-(7*Ffr.^2+8*Ffr+9)).*Alpha2 + ((7*Ffr.^2+8*Ffr+9)./((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1))).*atan(((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1)).*Alpha2);
% Ffr2 = (2*Ffr.^2+3*Ffr).*sin(1*atan(((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1)).*((1-(7*Ffr.^2+8*Ffr+9)).*Alpha2 + ((7*Ffr.^2+8*Ffr+9)./((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1))).*atan(((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1)).*Alpha2))))
%Ffr2 = (2*Ffr.^2+3*Ffr).*sin(1*atan(((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1)).*((1-(7*Ffr.^2+8*Ffr+9)).*Alpha2 + ((7*Ffr.^2+8*Ffr+9)./((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1))).*atan(((4*sin(5*atan(6*Ffr)))./((2*Ffr.^2+3*Ffr)*1)).*Alpha2))))


TireForce = @(data,slip,a)(a(2)*data.^2+a(3)*data).*sin(a(1)*atan(((a(4)*sin(a(5)*atan(a(6)*data)))./((a(2)*data.^2+a(3)*data)*a(1))).*((1-(a(7)*data.^2+a(8)*data+a(9))).*(slip) + ((a(7)*data.^2+a(8)*data+a(9))./((a(4)*sin(a(5)*atan(a(6)*data)))./((a(2)*data.^2+a(3)*data)*a(1)))).*atan(((a(4)*sin(a(5)*atan(a(6)*data)))./((a(2)*data.^2+a(3)*data)*a(1))).*(slip)))))
FTotFcn = @(a,DataX) TireForce(DataX(mI:mF,5),DataX(mI:mF,1),a)+TireForce(DataX(mI:mF,6),DataX(mI:mF,2),a)+TireForce(DataX(mI:mF,7),DataX(mI:mF,3),a)+TireForce(DataX(mI:mF,8),DataX(mI:mF,4),a)
                    
    
                  % + (a(2)*DataX(mI:mF,6).^2+a(3)*DataX(mI:mF,6)).*sin(a(1)*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,6))))./((a(2)*DataX(mI:mF,6).^2+a(3)*DataX(mI:mF,6))*a(1))).*((1-(a(7)*DataX(mI:mF,6).^2+a(8)*DataX(mI:mF,6)+a(9))).*(DataX(mI:mF,2)) + ((a(7)*DataX(mI:mF,6).^2+a(8)*DataX(mI:mF,6)+a(9))./((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,6))))./((a(2)*DataX(mI:mF,6).^2+a(3)*DataX(mI:mF,6))*a(1)))).*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,6))))./((a(2)*DataX(mI:mF,6).^2+a(3)*DataX(mI:mF,6))*a(1))).*(DataX(mI:mF,2))))))...
                  % + (a(2)*DataX(mI:mF,7).^2+a(3)*DataX(mI:mF,7)).*sin(a(1)*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,7))))./((a(2)*DataX(mI:mF,7).^2+a(3)*DataX(mI:mF,7))*a(1))).*((1-(a(7)*DataX(mI:mF,7).^2+a(8)*DataX(mI:mF,7)+a(9))).*(DataX(mI:mF,3)) + ((a(7)*DataX(mI:mF,7).^2+a(8)*DataX(mI:mF,7)+a(9))./((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,7))))./((a(2)*DataX(mI:mF,7).^2+a(3)*DataX(mI:mF,7))*a(1)))).*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,7))))./((a(2)*DataX(mI:mF,7).^2+a(3)*DataX(mI:mF,7))*a(1))).*(DataX(mI:mF,3))))))...
                  % + (a(2)*DataX(mI:mF,8).^2+a(3)*DataX(mI:mF,8)).*sin(a(1)*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,8))))./((a(2)*DataX(mI:mF,8).^2+a(3)*DataX(mI:mF,8))*a(1))).*((1-(a(7)*DataX(mI:mF,8).^2+a(8)*DataX(mI:mF,8)+a(9))).*(DataX(mI:mF,4)) + ((a(7)*DataX(mI:mF,8).^2+a(8)*DataX(mI:mF,8)+a(9))./((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,8))))./((a(2)*DataX(mI:mF,8).^2+a(3)*DataX(mI:mF,8))*a(1)))).*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,8))))./((a(2)*DataX(mI:mF,8).^2+a(3)*DataX(mI:mF,8))*a(1))).*(DataX(mI:mF,4))))));
%FTotFcn = @(a,DataX)(a(2)*DataX(mI:mF,5).^2+a(3)*DataX(mI:mF,5)).*sin(a(1)*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,5))))./((a(2)*DataX(mI:mF,5).^2+a(3)*DataX(mI:mF,5))*a(1))).*((1-(a(7)*DataX(mI:mF,5).^2+a(8)*DataX(mI:mF,5)+a(9))).*DataX(mI:mF,9) + ((a(7)*DataX(mI:mF,5).^2+a(8)*DataX(mI:mF,5)+a(9))./((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,5))))./((a(2)*DataX(mI:mF,5).^2+a(3)*DataX(mI:mF,5))*a(1)))).*atan(((a(4)*sin(a(5)*atan(a(6)*DataX(mI:mF,5))))./((a(2)*DataX(mI:mF,5).^2+a(3)*DataX(mI:mF,5))*a(1))).*DataX(mI:mF,9)))))

% C = a(1)
% D = Ffr*mu
% BCD = a(2)*sin(a(3)*atan(a(4)*Ffr))
% B = (a(2)*sin(a(3)*atan(a(4)*Ffr)))/(a(1)*Ffr*mu)
% E = a(5)*(Ffr)^2+a(6)*Ffr+a(7)


%Fy = (Ffr*.95).*sin((1)*atan(((2*sin(3*atan(4.*Ffr)))./(1*.95*Ffr)).*Alpha2-(5*(Ffr).^2+6*Ffr+7).*(((2*sin(3*atan(4.*Ffr)))./(1*.95*Ffr)).*Alpha2-atan(((2*sin(3*atan(4.*Ffr)))./(1*.95*Ffr)).*Alpha2))))
% Fy = (DataX(:,5)*.95).*sin((1)*atan(((2*sin(3*atan(4.*DataX(:,5))))./(1*.95*DataX(:,5))).*DataX(:,1)-(5*(DataX(:,5)).^2+6*DataX(:,5)+7).*(((2*sin(3*atan(4.*DataX(:,5))))./(1*.95*DataX(:,5))).*DataX(:,1)-atan(((2*sin(3*atan(4.*DataX(:,5))))./(1*.95*DataX(:,5))).*DataX(:,1)))))...
%    + (DataX(:,6)*.95).*sin((1)*atan(((2*sin(3*atan(4.*DataX(:,6))))./(1*.95*DataX(:,6))).*DataX(:,2)-(5*(DataX(:,6)).^2+6*DataX(:,6)+7).*(((2*sin(3*atan(4.*DataX(:,6))))./(1*.95*DataX(:,6))).*DataX(:,2)-atan(((2*sin(3*atan(4.*DataX(:,6))))./(1*.95*DataX(:,6))).*DataX(:,2)))))...
%    + (DataX(:,7)*.95).*sin((1)*atan(((2*sin(3*atan(4.*DataX(:,7))))./(1*.95*DataX(:,7))).*DataX(:,3)-(5*(DataX(:,7)).^2+6*DataX(:,7)+7).*(((2*sin(3*atan(4.*DataX(:,7))))./(1*.95*DataX(:,7))).*DataX(:,3)-atan(((2*sin(3*atan(4.*DataX(:,7))))./(1*.95*DataX(:,7))).*DataX(:,3)))))...
%    + (DataX(:,8)*.95).*sin((1)*atan(((2*sin(3*atan(4.*DataX(:,8))))./(1*.95*DataX(:,8))).*DataX(:,4)-(5*(DataX(:,8)).^2+6*DataX(:,8)+7).*(((2*sin(3*atan(4.*DataX(:,8))))./(1*.95*DataX(:,8))).*DataX(:,4)-atan(((2*sin(3*atan(4.*DataX(:,8))))./(1*.95*DataX(:,8))).*DataX(:,4)))))


%Fxfr = Ffr*mu*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4)*Ffr)))/(a(1)*Ffr*mu))*Alpha2-(a(5)*(Ffr)^2+a(6)*Ffr+a(7)*(((a(2)*sin(a(3)*atan(a(4)*Ffr)))/(a(1)*Ffr*mu))*Alpha2-atan(((a(2)*sin(a(3)*atan(a(4)*Ffr)))/(a(1)*Ffr*mu))*Alpha2)))))

% FTotFcn = @(a,Alpha2,Alpha3,Alpha4,Alpha5,Ffr,Ffl,Frr,Frl,mu)(Ffr*mu*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4)*Ffr)))/(a(1)*Ffr*mu))*Alpha2-(a(5)*(Ffr)^2+a(6)*Ffr+a(7)*(((a(2)*sin(a(3)*atan(a(4)*Ffr)))/(a(1)*Ffr*mu))*Alpha2-atan(((a(2)*sin(a(3)*atan(a(4)*Ffr)))/(a(1)*Ffr*mu))*Alpha2))))))...
%                                                            + (Ffl*mu*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4)*Ffl)))/(a(1)*Ffl*mu))*Alpha3-(a(5)*(Ffl)^2+a(6)*Ffl+a(7)*(((a(2)*sin(a(3)*atan(a(4)*Ffl)))/(a(1)*Ffl*mu))*Alpha3-atan(((a(2)*sin(a(3)*atan(a(4)*Ffl)))/(a(1)*Ffl*mu))*Alpha3))))))...
%                                                            + (Frr*mu*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4)*Frr)))/(a(1)*Frr*mu))*Alpha4-(a(5)*(Frr)^2+a(6)*Frr+a(7)*(((a(2)*sin(a(3)*atan(a(4)*Frr)))/(a(1)*Frr*mu))*Alpha4-atan(((a(2)*sin(a(3)*atan(a(4)*Frr)))/(a(1)*Frr*mu))*Alpha4))))))...
%                                                            + (Frl*mu*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4)*Frl)))/(a(1)*Frl*mu))*Alpha5-(a(5)*(Frl)^2+a(6)*Frl+a(7)*(((a(2)*sin(a(3)*atan(a(4)*Frl)))/(a(1)*Frl*mu))*Alpha5-atan(((a(2)*sin(a(3)*atan(a(4)*Frl)))/(a(1)*Frl*mu))*Alpha5))))));
% FTotFcn = @(a,DataX)(DataX(:,5).*.95*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,5))))/(a(1).*DataX(:,5).*.95)).*DataX(:,1)-(a(5).*(DataX(:,5)).^2+a(6)*DataX(:,5)+a(7)*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,5))))/(a(1).*DataX(:,5).*.95)).*DataX(:,1)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,5))))/(a(1).*DataX(:,5).*.95)).*DataX(:,1)))))))...
%                   + (DataX(:,6).*.95*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,6))))/(a(1).*DataX(:,6).*.95)).*DataX(:,2)-(a(5).*(DataX(:,6)).^2+a(6)*DataX(:,6)+a(7)*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,6))))/(a(1).*DataX(:,6).*.95)).*DataX(:,2)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,6))))/(a(1).*DataX(:,6).*.95)).*DataX(:,2)))))))...
%                   + (DataX(:,7).*.95*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,7))))/(a(1).*DataX(:,7).*.95)).*DataX(:,3)-(a(5).*(DataX(:,7)).^2+a(6)*DataX(:,7)+a(7)*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,7))))/(a(1).*DataX(:,7).*.95)).*DataX(:,3)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,7))))/(a(1).*DataX(:,7).*.95)).*DataX(:,3)))))))...
%                   + (DataX(:,8).*.95*sin(a(1)*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,8))))/(a(1).*DataX(:,8).*.95)).*DataX(:,4)-(a(5).*(DataX(:,8)).^2+a(6)*DataX(:,8)+a(7)*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,8))))/(a(1).*DataX(:,8).*.95)).*DataX(:,4)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,8))))/(a(1).*DataX(:,8).*.95)).*DataX(:,4)))))));
% 
% FTotFcn = @(a,DataX)(DataX(:,5)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,5))))./(a(1)*.95*DataX(:,5))).*DataX(:,1)-(a(5)*(DataX(:,5)).^2+a(6)*DataX(:,5)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,5))))./(a(1)*.95*DataX(:,5))).*DataX(:,1)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,5))))./(a(1)*.95*DataX(:,5))).*DataX(:,1)))))...
%                   + (DataX(:,6)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,6))))./(a(1)*.95*DataX(:,6))).*DataX(:,2)-(a(5)*(DataX(:,6)).^2+a(6)*DataX(:,6)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,6))))./(a(1)*.95*DataX(:,6))).*DataX(:,2)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,6))))./(a(1)*.95*DataX(:,6))).*DataX(:,2)))))...
%                   + (DataX(:,7)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,7))))./(a(1)*.95*DataX(:,7))).*DataX(:,3)-(a(5)*(DataX(:,7)).^2+a(6)*DataX(:,7)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,7))))./(a(1)*.95*DataX(:,7))).*DataX(:,3)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,7))))./(a(1)*.95*DataX(:,7))).*DataX(:,3)))))...
%                   + (DataX(:,8)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,8))))./(a(1)*.95*DataX(:,8))).*DataX(:,4)-(a(5)*(DataX(:,8)).^2+a(6)*DataX(:,8)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(:,8))))./(a(1)*.95*DataX(:,8))).*DataX(:,4)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(:,8))))./(a(1)*.95*DataX(:,8))).*DataX(:,4)))));

% FTotFcn = @(a,DataX)(DataX(6001:15701,5)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,5))))./(a(1)*.95*DataX(6001:15701,5))).*DataX(6001:15701,1)-(a(5)*(DataX(6001:15701,5)).^2+a(6)*DataX(6001:15701,5)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,5))))./(a(1)*.95*DataX(6001:15701,5))).*DataX(6001:15701,1)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,5))))./(a(1)*.95*DataX(6001:15701,5))).*DataX(6001:15701,1)))))...
%                   + (DataX(6001:15701,6)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,6))))./(a(1)*.95*DataX(6001:15701,6))).*DataX(6001:15701,2)-(a(5)*(DataX(6001:15701,6)).^2+a(6)*DataX(6001:15701,6)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,6))))./(a(1)*.95*DataX(6001:15701,6))).*DataX(6001:15701,2)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,6))))./(a(1)*.95*DataX(6001:15701,6))).*DataX(6001:15701,2)))))...
%                   + (DataX(6001:15701,7)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,7))))./(a(1)*.95*DataX(6001:15701,7))).*DataX(6001:15701,3)-(a(5)*(DataX(6001:15701,7)).^2+a(6)*DataX(6001:15701,7)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,7))))./(a(1)*.95*DataX(6001:15701,7))).*DataX(6001:15701,3)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,7))))./(a(1)*.95*DataX(6001:15701,7))).*DataX(6001:15701,3)))))...
%                   + (DataX(6001:15701,8)*.95).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,8))))./(a(1)*.95*DataX(6001:15701,8))).*DataX(6001:15701,4)-(a(5)*(DataX(6001:15701,8)).^2+a(6)*DataX(6001:15701,8)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,8))))./(a(1)*.95*DataX(6001:15701,8))).*DataX(6001:15701,4)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,8))))./(a(1)*.95*DataX(6001:15701,8))).*DataX(6001:15701,4)))));
% FTotFcn = @(a,DataX)(a(8)*(DataX(6001:15701,5)).^2+a(9)*(DataX(6001:15701,5))).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,5))))./(a(1)*.95*DataX(6001:15701,5))).*DataX(6001:15701,1)-(a(5)*(DataX(6001:15701,5)).^2+a(6)*DataX(6001:15701,5)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,5))))./(a(1)*.95*DataX(6001:15701,5))).*DataX(6001:15701,1)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,5))))./(a(1)*.95*DataX(6001:15701,5))).*DataX(6001:15701,1)))))...
%                   + (a(8)*(DataX(6001:15701,6)).^2+a(9)*(DataX(6001:15701,6))).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,6))))./(a(1)*.95*DataX(6001:15701,6))).*DataX(6001:15701,2)-(a(5)*(DataX(6001:15701,6)).^2+a(6)*DataX(6001:15701,6)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,6))))./(a(1)*.95*DataX(6001:15701,6))).*DataX(6001:15701,2)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,6))))./(a(1)*.95*DataX(6001:15701,6))).*DataX(6001:15701,2)))))...
%                   + (a(8)*(DataX(6001:15701,7)).^2+a(9)*(DataX(6001:15701,7))).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,7))))./(a(1)*.95*DataX(6001:15701,7))).*DataX(6001:15701,3)-(a(5)*(DataX(6001:15701,7)).^2+a(6)*DataX(6001:15701,7)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,7))))./(a(1)*.95*DataX(6001:15701,7))).*DataX(6001:15701,3)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,7))))./(a(1)*.95*DataX(6001:15701,7))).*DataX(6001:15701,3)))))...
%                   + (a(8)*(DataX(6001:15701,8)).^2+a(9)*(DataX(6001:15701,8))).*sin((a(1))*atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,8))))./(a(1)*.95*DataX(6001:15701,8))).*DataX(6001:15701,4)-(a(5)*(DataX(6001:15701,8)).^2+a(6)*DataX(6001:15701,8)+a(7)).*(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,8))))./(a(1)*.95*DataX(6001:15701,8))).*DataX(6001:15701,4)-atan(((a(2)*sin(a(3)*atan(a(4).*DataX(6001:15701,8))))./(a(1)*.95*DataX(6001:15701,8))).*DataX(6001:15701,4)))));
% 

%a0 = [-.29,-.0743,2.2965,10.7338,9.7358,-.0814,4.5641,1.281,-5.8646]
% a0 = [.2483,.006,1.4468,5.2278,1.0563,.3656,-1.7396,.0436,5.983]
%a = [.45,-.08,2,8,1.0579,.3673,-5.3059,1.000257,5.983]
% a = [.3729,-.0364,1.8348,7.0561,1.0579,.3673,-3.3059,.00257,5.983]
%a0 = [.2701,-.0064,.9009,10,1.0571,.3672,-5.306,.0026,5.9830]
%a0 = [.3541,-.4076,.9521,61.9389,1.0193,1.0792,-3.9750,-.0522,5.9777]
%a0 = [.5258,-.0480,.9560,708.1228,1.9591,25.4768,-31.3957,.1,.1]
%a0 = [.0941,-.1845,.8396,10.2904,1.0574,.3679,-5.4828,-.0030,5.9828]
%a0 = [.2266,.1521,-6.2750,15,10,.0028,-40,-10,-10]
%a0 = [.8,.0005,1.1,70,1.85,385.01,-48.85,10,1]
%a0 = [1.1,-.001,1.1,70,.9994,299.9632,-5.2699,10,10]
%a0 = [1.1487,-.0096,1.1,110,49,-388,0.0127,33.65,10]
%a0 = [1.2508,.0010,1.05,52.3878,.092,5.0831,-16.6620,-10,-10]
%a0 = [1.1277,-.0015,1.1685,72.4772,.0678,68.2138,-10.2534,-161.4905,-15.8796]
%a0 = [.1,.1,.1,.1,.1,.1,.1,.1,.1]
a0 = [.8 .001 1.05 155 .6 6 -.002 0 .707]
A = [0 2500 50 0 0 0 0 0 0];
B = 100;
LB = [-100,-.002,-100,-100,-100,-100,-400,-10,-10];
UB = [1.5,.001,1.1,110,100,300,100,10,10];
%LB = [0,-.001,0,0,0,0,-30,0,0]
options = optimoptions('lsqcurvefit','MaxIterations',100,'MaxFunctionEvaluations',100,'FiniteDifferenceStepSize',.0001)
[a,residual] = lsqcurvefit(FTotFcn,a0,DataX,(TotalForce(mI:mF)),[],[],[],[],[],[],[],options)
FcnOut = FTotFcn(a0,DataX)

%plot(residual)
% plot(AlphaTime(41701:149544),[TotalForce(41701:149544),CurveFitForce(41701:149544)])
% title('Steady-State Cornering Test Curve Fit')
% xlabel('Time (s)')
% ylabel('Lateral Force (N)')
% legend('From Tire Model','Experimental')
% %ylim([0 150])
% xlim([40 150])


% a0 = [1,1,1,1,1,1,1]
% options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',100000000)
% a = lsqcurvefit(FTotFcn,a0,DataX,TotalForce,[],[],[],[],[],[],[],options,[]);

% figure;
% plot(XdS.Time, XdS.Variables);
% hold on
% plot(YdS.Time, YdS.Variables);
% hold off
%plot(PsiS)





%AccelData = smoothdata(table2array(InputData(:,3)), "movmean", 10);
%TimeData = table2array(InputData(:,1));
%SteerTS = timetable(SteerData, 'TimeStep', seconds(.01166))
% VelTS = timetable(VelData, 'TimeStep', seconds(.0114))
% AngTS = timetable(AngData, 'TimeStep', seconds(.01))
% AccelTS = timetable(AccelData, 'TimeStep', seconds(.01))
%SteerTS = timeseries(SteerData, TimeData);
%SteerTS2 = timetable(TimeData, SteerData);


