function [ DF_out ] = LG_DiffFlat(traj)
%% function to generate states and control from flat outputs
%% parameters:
mP = 0.2;
mQ1 = 0.55;
mQ2 = 0.55;
IA = diag([0.0023,0.0028,0.0046]);
IB = diag([0.0023,0.0028,0.0046]);
g = 9.8;
e3 = [0,0,1];
Lr = 1;
N = length(traj.time);

for i=1:N
    xP0 = traj.x(1,i);
    xP1 = traj.x(2,i);
    xP2 = traj.x(3,i);
    xP3 = traj.x(4,i);
    xP4 = traj.x(5,i);
    xP5 = traj.x(6,i);
    xP6 = traj.x(7,i);
    
    yP0 = traj.y(1,i);
    yP1 = traj.y(2,i);
    yP2 = traj.y(3,i);
    yP3 = traj.y(4,i);
    yP4 = traj.y(5,i);
    yP5 = traj.y(6,i);
    yP6 = traj.y(7,i);

    zP0 = traj.z(1,i);
    zP1 = traj.z(2,i);
    zP2 = traj.z(3,i);
    zP3 = traj.z(4,i);
    zP4 = traj.z(5,i);
    zP5 = traj.z(6,i);
    zP6 = traj.z(7,i);
    
    aB0 = traj.aB(1,i)+deg2rad(180);
    aB1 = traj.aB(2,i);    
    aB2 = traj.aB(3,i);    
    aB3 = traj.aB(4,i);
    aB4 = traj.aB(5,i);    
    aB5 = traj.aB(6,i);   
    aB6 = traj.aB(7,i); 
    
    DF_out.aA0(i) = traj.aB(1,i);              % alpha1
    
    
    DF_out.aB0(i) = traj.aB(1,i)+deg2rad(180); % alpha2
    DF_out.aB1(i) = traj.aB(2,i);              % dalpha2
    DF_out.bB0(i) = deg2rad(45);               % beta2
    DF_out.bB1(i) = 0;                         % dbeta2    
   %% payload position
    DF_out.xP(i) = xP0;
    DF_out.yP(i) = yP0;
    DF_out.zP(i) = zP0;
    DF_out.dxP(i) = xP1;
    DF_out.dyP(i) = yP1;
    DF_out.dzP(i) = zP1;
    DF_out.axP(i) = xP2;
    DF_out.ayP(i) = yP2;
    DF_out.azP(i) = zP2;
    TB = g*mP/sqrt(2);
    %% Tension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% roB,
    DF_out.roB0(:,i) = [2.^(-1/2).*cos(aB0),2.^(-1/2).*sin(aB0),2.^(-1/2)];
    DF_out.roB1(:,i) = [(-1).*2.^(-1/2).*aB1.*sin(aB0),2.^(-1/2).*aB1.*cos(aB0),0];
    DF_out.roB2(:,i) = [(-1).*2.^(-1/2).*(aB1.^2.*cos(aB0)+aB2.*sin(aB0)),2.^(-1/2).*( ...
  aB2.*cos(aB0)+(-1).*aB1.^2.*sin(aB0)),0];
    DF_out.roB3(:,i) = [2.^(-1/2).*((-3).*aB1.*aB2.*cos(aB0)+(aB1.^3+(-1).*aB3).*sin(aB0) ...
  ),2.^(-1/2).*(((-1).*aB1.^3+aB3).*cos(aB0)+(-3).*aB1.*aB2.*sin( ...
  aB0)),0];
    DF_out.roB4(:,i) = [2.^(-1/2).*((aB1.^4+(-3).*aB2.^2+(-4).*aB1.*aB3).*cos(aB0)+(6.* ...
  aB1.^2.*aB2+(-1).*aB4).*sin(aB0)),2.^(-1/2).*(((-6).*aB1.^2.*aB2+ ...
  aB4).*cos(aB0)+(aB1.^4+(-3).*aB2.^2+(-4).*aB1.*aB3).*sin(aB0)),0];
    DF_out.TBroB0(:,i) = TB*DF_out.roB0(:,i);
    DF_out.TBroB1(:,i) = TB*DF_out.roB1(:,i);
    DF_out.TBroB2(:,i) = TB*DF_out.roB2(:,i);
    DF_out.TBroB3(:,i) = TB*DF_out.roB3(:,i);
    DF_out.TBroB4(:,i) = TB*DF_out.roB4(:,i);
    

%     DF_out.
    
    
    
   %%% roA0,TA0,  aA0,bA0
   T1ro10 = mP*[xP2,yP2,zP2]'-DF_out.TBroB0(:,i)+mP*g*[0,0,1]';
   
    DF_out.roA0(:,i) = T1ro10/norm(T1ro10,2);
    %% alphaA,betaA
    [tempa,tempb] = LG_ro2ab(DF_out.roA0(:,i));
%     DF_out.aA0(i) = traj.aB(1,i);
    DF_out.aA0(i) = tempa;
    DF_out.bA0(i) = tempb;
    
    
    DF_out.TA0(i) = dot(T1ro10,DF_out.roA0(:,i));
    %%% roA1,TA1
    DF_out.TA1(i) = dot(mP*[xP3,yP3,zP3]',DF_out.roA0(:,i))-...
       dot(DF_out.TBroB1(:,i),DF_out.roA0(:,i));
    DF_out.roA1(:,i) = DF_out.TA0(i)\(mP*[xP3,yP3,zP3]'-DF_out.TA1(i)*DF_out.roA0(:,i)-DF_out.TBroB1(:,i));
    
    
    
    %% dalphaA,dbetaA
    roA1 = DF_out.roA1(:,i);
    DF_out.bA1(i) = roA1(3)/cos(DF_out.bA0(i));
    a1 = DF_out.aA0(i);
    b1 = DF_out.bA0(i);
    DF_out.aA1(i) = (roA1(1)+sin(b1)*cos(a1)*DF_out.bA1(i))/(-cos(b1)*sin(a1));
    
    %%% roA2,TA2
    DF_out.TA2(i) = dot(mP*[xP4,yP4,zP4]',DF_out.roA0(:,i)) +...
        dot(mP*[xP3,yP3,zP3]',DF_out.roA1(:,i)) - dot(DF_out.TBroB2(:,i),DF_out.roA0(:,i))-...
        dot(DF_out.TBroB1(:,i),DF_out.roA1(:,i));
    
    temp1 = mP*[xP4,yP4,zP4]';
    temp2 = DF_out.TA2(i)*DF_out.roA0(:,i);
    temp3 = 2*DF_out.TA1(i)*DF_out.roA1(:,i);
    DF_out.roA2(:,i) = DF_out.TA0(i)\(temp1-temp2-temp3-DF_out.TBroB2(:,i));
    %%% roA3,TA3
    temp1 = dot(mP*[xP5,yP5,zP5]',DF_out.roA0(:,i))+dot(mP*[xP4,yP4,zP4]',DF_out.roA1(:,i));
    temp2 = dot(mP*[xP4,yP4,zP4]',DF_out.roA1(:,i))+dot(mP*[xP3,yP3,zP3]',DF_out.roA2(:,i));
    temp3 = dot(DF_out.TBroB3(:,i),DF_out.roA0(:,i))+dot(DF_out.TBroB2(:,i),DF_out.roA1(:,i));
    temp4 = dot(DF_out.TBroB2(:,i),DF_out.roA1(:,i))+dot(DF_out.TBroB1(:,i),DF_out.roA2(:,i));
    DF_out.TA3(i) = (temp1+temp2-temp3-temp4);
    
    temp1 = mP*[xP5,yP5,zP5]';
    temp2 = DF_out.TA3(i)*DF_out.roA0(:,i);
    temp3 = 3*DF_out.TA1(i)*DF_out.roA2(:,i);
    temp4 = 3*DF_out.TA2(i)*DF_out.roA1(:,i);
    temp5 = DF_out.TBroB3(:,i);
    DF_out.roA3(:,i) = DF_out.TA0(i)\(temp1-temp2-temp3-temp4-temp5);
    
    %%%  roA4,TA4
    temp1 = dot(mP*[xP6,yP6,zP6]',DF_out.roA0(:,i));
    temp2 = 3*dot(mP*[xP5,yP5,zP5]',DF_out.roA1(:,i));
    temp3 = 3*dot(mP*[xP4,yP4,zP4]',DF_out.roA2(:,i));
    temp4 = dot(mP*[xP3,yP3,zP3]',DF_out.roA3(:,i));
    temp5 = dot(DF_out.TBroB4(:,i),DF_out.roA0(:,i));
    temp6 = 3*dot(DF_out.TBroB3(:,i),DF_out.roA1(:,i));
    temp7 = 3*dot(DF_out.TBroB2(:,i),DF_out.roA2(:,i));
    temp8 = dot(DF_out.TBroB1(:,i),DF_out.roA3(:,i));
    
    DF_out.TA4(i) = temp1+temp2+temp3+temp4-temp5-temp6-temp7-temp8;
    
    temp1 = mP*[xP6,yP6,zP6]';
    temp2 =   DF_out.TA4(i)*DF_out.roA0(:,i);
    temp3 = 4*DF_out.TA3(i)*DF_out.roA1(:,i);
    temp4 = 6*DF_out.TA2(i)*DF_out.roA2(:,i);
    temp5 = 4*DF_out.TA1(i)*DF_out.roA3(:,i);
    temp6 = DF_out.TBroB4(:,i);
    
    DF_out.roA4(:,i) =DF_out.TA0(i)\(temp1-temp2-temp3-temp4-temp5-temp6);
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% quad1 position
    tempPosQ1 = [xP0;yP0;zP0] + Lr*DF_out.roA0(:,i);
    
    DF_out.xQ1(i) = tempPosQ1(1);
    DF_out.yQ1(i) = tempPosQ1(2);
    DF_out.zQ1(i) = tempPosQ1(3);

    tempVelQ1 = [xP1;yP1;zP1] + Lr*DF_out.roA1(:,i);
    DF_out.dxQ1(i) = tempVelQ1(1);
    DF_out.dyQ1(i) = tempVelQ1(2);
    DF_out.dzQ1(i) = tempVelQ1(3);

    tempAccQ1 = [xP2;yP2;zP2] + Lr*DF_out.roA2(:,i);
    
    DF_out.axQ1(i) = tempAccQ1(1);
    DF_out.ayQ1(i) = tempAccQ1(2);
    DF_out.azQ1(i) = tempAccQ1(3);
    
    tempJekQ1 = [xP3;yP3;zP3] + Lr*DF_out.roA3(:,i);
    
    DF_out.jxQ1(i) = tempJekQ1(1);
    DF_out.jyQ1(i) = tempJekQ1(2);
    DF_out.jzQ1(i) = tempJekQ1(3);
    
    tempSnpQ1 = [xP4;yP4;zP4] + Lr*DF_out.roA4(:,i);
    
    DF_out.sxQ1(i) = tempSnpQ1(1);
    DF_out.syQ1(i) = tempSnpQ1(2);
    DF_out.szQ1(i) = tempSnpQ1(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% quad2 position
%     tempPosQ2 = [xP0+2.^(-1/2).*Lr.*cos(aB0),yP0+2.^(-1/2).*Lr.*sin(aB0),2.^(-1/2) ...
%   .*Lr+zP0];
    tempPosQ2 = [xP0;yP0;zP0] + Lr*DF_out.roB0(:,i);
    
    DF_out.xQ2(i) = tempPosQ2(1);
    DF_out.yQ2(i) = tempPosQ2(2);
    DF_out.zQ2(i) = tempPosQ2(3);
    
    tempVelQ2 = [xP1;yP1;zP1] + Lr*DF_out.roB1(:,i);
    %[xP1+(-1).*2.^(-1/2).*aB1.*Lr.*sin(aB0),yP1+2.^(-1/2).*aB1.*Lr.* ...
%   cos(aB0),zP1];

    DF_out.dxQ2(i) = tempVelQ2(1);
    DF_out.dyQ2(i) = tempVelQ2(2);
    DF_out.dzQ2(i) = tempVelQ2(3);
    
    tempAccQ2 = [xP2;yP2;zP2] + Lr*DF_out.roB2(:,i);
    %[xP2+(-1).*2.^(-1/2).*aB1.^2.*Lr.*cos(aB0)+(-1).*2.^(-1/2).*aB2.* ...
%   Lr.*sin(aB0),yP2+2.^(-1/2).*aB2.*Lr.*cos(aB0)+(-1).*2.^(-1/2).* ...
%   aB1.^2.*Lr.*sin(aB0),zP2]; 

    DF_out.axQ2(i) = tempAccQ2(1);
    DF_out.ayQ2(i) = tempAccQ2(2);
    DF_out.azQ2(i) = tempAccQ2(3);
    
    tempJekQ2 = [xP3;yP3;zP3] + Lr*DF_out.roB3(:,i);
    
    DF_out.jxQ2(i) = tempJekQ2(1);
    DF_out.jyQ2(i) = tempJekQ2(2);
    DF_out.jzQ2(i) = tempJekQ2(3);
    
    tempSnpQ2 = [xP4;yP4;zP4] + Lr*DF_out.roB4(:,i);
    
    DF_out.sxQ2(i) = tempSnpQ2(1);
    DF_out.syQ2(i) = tempSnpQ2(2);
    DF_out.szQ2(i) = tempSnpQ2(3);    
    
    
%% RA, fA, RB, fB
  f1zbA = mQ1*[DF_out.axQ1(i);DF_out.ayQ1(i);DF_out.azQ1(i)]+...
      mQ1*g*[0;0;1] + DF_out.TA0(i)*DF_out.roA0(:,i);
  zbA = f1zbA/norm(f1zbA,2);
  DF_out.fA(i) = dot(f1zbA,zbA);
  xcA = [cos(0);sin(0);0];
  ybA = cross(zbA,xcA)/norm(cross(zbA,xcA),2);
  xbA = cross(ybA,zbA);
  RA = [xbA,ybA,zbA];
  DF_out.eulerA(:,i) = rotMat2euler(RA);
  
  f2zbB = mQ2*[DF_out.axQ2(i);DF_out.ayQ2(i);DF_out.azQ2(i)]+...
      mQ2*g*[0;0;1] + DF_out.TBroB0(:,i);
  zbB = f2zbB/norm(f2zbB,2);
  DF_out.fB(i) = dot(f2zbB,zbB);
  xcB = [cos(0);sin(0);0];
  ybB = cross(zbB,xcB)/norm(cross(zbB,xcB),2);
  xbB = cross(ybB,zbB);
  RB = [xbB,ybB,zbB];
  DF_out.eulerB(:,i) = rotMat2euler(RB);
  
  %%%%%%%%%%%%%%%%%%%%%%%
  %% w1
  jerkQ1 = [DF_out.jxQ1(i);DF_out.jyQ1(i);DF_out.jzQ1(i)];
  dT1ro1 = DF_out.TA1(i)*DF_out.roA0(:,i)+DF_out.TA0(i)*DF_out.roA1(:,i);
  df1 = mQ1*dot(jerkQ1,zbA) + dot(dT1ro1,zbA);
  hw1 = (mQ1*jerkQ1 - df1*zbA + dT1ro1)/DF_out.fA(i);
  p1 = -dot(hw1,ybA);
  q1 = dot(hw1,xbA);
  r1 = 0;
  DF_out.w1(:,i)=[p1,q1,r1]';
%   DF_out.w1(:,i)= 
  %% w2
  jerkQ2 = [DF_out.jxQ2(i);DF_out.jyQ2(i);DF_out.jzQ2(i)];
  dT2ro2 = DF_out.TBroB1(:,i);
  df2 = mQ2*dot(jerkQ2,zbB) + dot(dT1ro1,zbB);
  hw2 = (mQ2*jerkQ2 - df2*zbB + dT2ro2)/DF_out.fB(i);
  p2 = -dot(hw2,ybB);
  q2 = dot(hw2,xbB);
  r2 = 0;
  DF_out.w2(:,i)=[p2,q2,r2]';%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %% dw1
  snapQ1 = [DF_out.sxQ1(i);DF_out.syQ1(i);DF_out.szQ1(i)];
  d2T1ro1 = DF_out.TA2(i)*DF_out.roA0(:,i) +...
      2*DF_out.TA1(i)*DF_out.roA1(:,i) + DF_out.TA0(i)*DF_out.roA2(:,i);
  temp1 = cross(DF_out.w1(:,i),DF_out.fA(i)*zbA);
  temp2 = cross(DF_out.w1(:,i),temp1);
  temp = dot(temp2,zbA);
  d2f1 = mQ1*dot(snapQ1,zbA)-temp + dot(d2T1ro1,zbA);
  
  temp3 = d2f1*zbA + cross(DF_out.w1(:,i),df1*zbA);
  temp4 = cross(DF_out.w1(:,i),df1*zbA+cross(DF_out.w1(:,i),DF_out.fA(i)*zbA));
  
  hdw1 = (mQ1*snapQ1 - temp3 - temp4 + d2T1ro1)/DF_out.fA(i);
  dp1 = -dot(hdw1,ybA);
  dq1 = dot(hdw1,xbA);
  dr1 = 0;
  DF_out.dw1(:,i)=[dp1,dq1,dr1]';
  %%%%%%%%%%%%%%%%%%%%%%%
  %% dw2
  snapQ2 = [DF_out.sxQ2(i);DF_out.syQ2(i);DF_out.szQ2(i)];
  d2T2ro2 = DF_out.TBroB2(:,i);
  temp1 = cross(DF_out.w2(:,i),DF_out.fB(i)*zbB);
  temp2 = cross(DF_out.w2(:,i),temp1);
  temp = dot(temp2,zbB);
  d2f2 = mQ2*dot(snapQ2,zbB)-temp + dot(d2T2ro2,zbB);
  
  temp3 = d2f2*zbB + cross(DF_out.w2(:,i),df2*zbB);
  temp4 = cross(DF_out.w2(:,i),df2*zbB+cross(DF_out.w2(:,i),DF_out.fB(i)*zbB));
  
  hdw2 = (mQ2*snapQ2 - temp3 - temp4 + d2T2ro2)/DF_out.fB(i);
  dp2 = -dot(hdw2,ybB);
  dq2 = dot(hdw2,xbB);
  dr2 = 0;
  
  DF_out.dw2(:,i)=[dp2,dq2,dr2]';  %%
  %%%
  DF_out.MA(:,i) = cross(DF_out.w1(:,i),IA*DF_out.w1(:,i)) + IA*DF_out.dw1(:,i);
  DF_out.MB(:,i) = cross(DF_out.w2(:,i),IB*DF_out.w2(:,i)) + IB*DF_out.dw2(:,i);
  
end

test = 1;
end
