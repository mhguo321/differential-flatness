function fig = LG_6DoFAnimation_door( varargin )
%LG_6DOFANIMATION Summary of this function goes here
%   Detailed explanation goes here
%% Inputs:
if(nargin==10)
    pQ1_total = varargin{1};
    R1_total = varargin{2};
    pQ2_total = varargin{3};
    R2_total  = varargin{4};
    pP_total = varargin{5};
    SamplePlotFreq = varargin{6};
    Trail = varargin{7};
    CreateAVI = varargin{8};
    isFixView = varargin{9};
    DF_out = varargin{10};
else
    pQ1_total = varargin{1};
    R1_total = varargin{2};
    pQ2_total = varargin{3};
    R2_total  = varargin{4};
    pP_total = varargin{5};
    SamplePlotFreq = varargin{6};
    Trail = varargin{7};
    CreateAVI = varargin{8};
    isFixView = varargin{9};   
end
[numSamples dummy] = size(pQ1_total);% total number of samples
%% Default configuration parameters:
FullScreen = false;
% Trail = 'DotsOnly';
AVIfileName = 'Demo';
AVIfileNameEnum = true; 

ShowArrowHead = true;
    if(ShowArrowHead)
        ShowArrowHeadStr = 'on';
    else
        ShowArrowHeadStr = 'off';
    end
    
ShowLegend = false;
Title = 'Quadcopters carrying a payload--';
Position = [50 50 866 600];
Spin = 120;
% View = [(-0:(Spin/(length(pP_total)-1)):(-0+Spin))', 0*ones(length(pP_total), 1)];
View = [(18:(Spin/(length(pP_total)-1)):(18+Spin))', 35*ones(length(pP_total), 1)];
% View = [(-0:(Spin/(length(pP_total)-1)):(-0+Spin))', 90*ones(length(pP_total), 1)];
AxisLength = 1.3;
LimitRatio = 1;
Xlabel = '$X(m)$';
Ylabel = '$Y(m)$';
Zlabel = '$Z(m)$';

pQ1 = pQ1_total(1:SamplePlotFreq:numSamples, :);
pQ2 = pQ2_total(1:SamplePlotFreq:numSamples, :);
R1 = R1_total(:, :, 1:SamplePlotFreq:numSamples) * AxisLength;
R2 = R2_total(:, :, 1:SamplePlotFreq:numSamples) * AxisLength;
pP = pP_total(1:SamplePlotFreq:numSamples, :);
if(numel(View) > 2)
    View = View(1:SamplePlotFreq:numSamples, :);
end
[numPlotSamples dummy] = size(pP);
    
%% Setup AVI file

    aviobj = [];                                                            	% create null object
    if(CreateAVI)
        fileName = strcat(AVIfileName, '.avi');
        if(exist(fileName, 'file'))
            if(AVIfileNameEnum)                                              	% if file name exists and enum enabled
                i = 0;
                while(exist(fileName, 'file'))                                  % find un-used file name by appending enum
                    fileName = strcat(AVIfileName, sprintf('%i', i), '.avi');
                    i = i + 1;
                end
            else                                                                % else file name exists and enum disabled
                fileName = [];                                                  % file will not be created
            end
        end
        if(isempty(fileName))
            sprintf('AVI file not created as file already exists.')
        else
%             aviobj = avifile(fileName, 'fps', AVIfps, 'compression', 'Cinepak', 'quality', 100,'VideoCompressionMethod','H.264');
              aviobj = VideoWriter(fileName);
              open(aviobj);
        end
    end    
%%
fig = figure( 'Name', '6DOF Animation(Email:gmh_njust@163.com)');
    if(FullScreen)
        screenSize = get(0, 'ScreenSize');
        set(fig, 'Position', [0 0 screenSize(3) screenSize(4)]);
    elseif(~isempty(Position))
        set(fig, 'Position', Position);
    end
    lighting phong;
    set(gcf, 'Renderer', 'zbuffer');
    hold on;
    axis equal;
    grid on;
    view([View(1,1),View(1,2)]);
%     title(i);
    xlabel(Xlabel,'Interpreter','Latex');
    ylabel(Ylabel,'Interpreter','Latex');
    zlabel(Zlabel,'Interpreter','Latex');
    
   % Create plot data arrays
    if(strcmp(Trail, 'DotsOnly') || strcmp(Trail, 'All'))
        xQ1 = zeros(numPlotSamples, 1);
        yQ1 = zeros(numPlotSamples, 1);
        zQ1 = zeros(numPlotSamples, 1);
        
        xQ2 = zeros(numPlotSamples, 1);
        yQ2 = zeros(numPlotSamples, 1);
        zQ2 = zeros(numPlotSamples, 1);
        
        xP = zeros(numPlotSamples, 1);
        yP = zeros(numPlotSamples, 1);
        zP = zeros(numPlotSamples, 1); 
    end
    if(strcmp(Trail, 'All'))
        xQ1 = zeros(numPlotSamples, 1);
        yQ1 = zeros(numPlotSamples, 1);
        zQ1 = zeros(numPlotSamples, 1);
        
        xQ2 = zeros(numPlotSamples, 1);
        yQ2 = zeros(numPlotSamples, 1);
        zQ2 = zeros(numPlotSamples, 1);
        
        xP = zeros(numPlotSamples, 1);
        yP = zeros(numPlotSamples, 1);
        zP = zeros(numPlotSamples, 1);
        
        ux1 = zeros(numPlotSamples, 1);
        vx1 = zeros(numPlotSamples, 1);
        wx1 = zeros(numPlotSamples, 1);
        uy1 = zeros(numPlotSamples, 1);
        vy1 = zeros(numPlotSamples, 1);
        wy1 = zeros(numPlotSamples, 1);
        uz1 = zeros(numPlotSamples, 1);
        vz1 = zeros(numPlotSamples, 1);
        wz1 = zeros(numPlotSamples, 1);
        
        ux2 = zeros(numPlotSamples, 1);
        vx2 = zeros(numPlotSamples, 1);
        wx2 = zeros(numPlotSamples, 1);
        uy2 = zeros(numPlotSamples, 1);
        vy2 = zeros(numPlotSamples, 1);
        wy2 = zeros(numPlotSamples, 1);
        uz2 = zeros(numPlotSamples, 1);
        vz2 = zeros(numPlotSamples, 1);
        wz2 = zeros(numPlotSamples, 1);
    end
    xQ1(1) = pQ1(1,1);
    yQ1(1) = pQ1(1,2);
    zQ1(1) = pQ1(1,3);  
    
    xQ2(1) = pQ2(1,1);
    yQ2(1) = pQ2(1,2);
    zQ2(1) = pQ2(1,3);  
    
    xP(1) = pP(1,1);
    yP(1) = pP(1,2);
    zP(1) = pP(1,3);  
    
    oxQ1(1) = pQ1(1,1);% Quadcopter1 initial position
    oyQ1(1) = pQ1(1,2);
    ozQ1(1) = pQ1(1,3);
    
    oxQ2(1) = pQ2(1,1);% Quadcopter2 initial position
    oyQ2(1) = pQ2(1,2);
    ozQ2(1) = pQ2(1,3);
    
    oxP(1) = pP(1,1);% payload initial position
    oyP(1) = pP(1,2);
    ozP(1) = pP(1,3);
    
  
%     ratio_quiver = 1;
%     R = ratio_quiver.*R;
    ux1(1) = R1(1,1,1:1);% each column vector of R is a direction vector
    vx1(1) = R1(2,1,1:1);
    wx1(1) = R1(3,1,1:1);
    uy1(1) = R1(1,2,1:1);
    vy1(1) = R1(2,2,1:1);
    wy1(1) = R1(3,2,1:1);
    uz1(1) = R1(1,3,1:1);
    vz1(1) = R1(2,3,1:1);
    wz1(1) = R1(3,3,1:1); 
    
    ux2(1) = R2(1,1,1:1);% each column vector of R is a direction vector
    vx2(1) = R2(2,1,1:1);
    wx2(1) = R2(3,1,1:1);
    uy2(1) = R2(1,2,1:1);
    vy2(1) = R2(2,2,1:1);
    wy2(1) = R2(3,2,1:1);
    uz2(1) = R2(1,3,1:1);
    vz2(1) = R2(2,3,1:1);
    wz2(1) = R2(3,3,1:1); 


    
    % Create graphics handles
%     subplot(2,1,1);
    orginalHandle = plot3(0,0,0,'rh');
    if(nargin==10)
        trajHandle = plot3(DF_out.xP,DF_out.yP,DF_out.zP,'r-.');
    end
    line1XLhandle = line('XData',[2-0.5,2-0.5],'YData',[3-0.0,3-0.0],'ZData',0.8+[-0.2,1.3],'LineWidth',4,'Color',[128 42 42]/256.0);

    line1ZDhandle = line('XData',[2-0.5,2+0.5],'YData',[3-0.0,3+0.0],'ZData',0.8+[-0.2,-0.2],'LineWidth',4,'Color',[128 42 42]/256.0);
    
    org1Handle = plot3(oxQ1, oyQ1, ozQ1, 'k.'); % mass center of the quadcopter1
    
    
    quivX1handle = quiver3(oxQ1, oyQ1, ozQ1, ux1, vx1, wx1,  'r', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    quivY1handle = quiver3(oxQ1, oyQ1, ozQ1, uy1, vy1, wy1,  'g', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    quivZ1handle = quiver3(oxQ1, oyQ1, ozQ1, uz1, vz1, wz1,  'b', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    
    len1 = 0.34; % distance from the center of rotor 1 to rotor 3; 
    Q1x1 = [oxQ1,oyQ1,ozQ1] + 0.5*len1*[ux1,vx1,wx1];
    Q1x3 = [oxQ1,oyQ1,ozQ1] - 0.5*len1*[ux1,vx1,wx1];
    Q1x2 = [oxQ1,oyQ1,ozQ1] + 0.5*len1*[uy1,vy1,wy1];
    Q1x4 = [oxQ1,oyQ1,ozQ1] - 0.5*len1*[uy1,vy1,wy1];
    Q1X1 = [Q1x1(1),Q1x3(1)];
    Q1X2 = [Q1x1(2),Q1x3(2)];
    Q1X3 = [Q1x1(3),Q1x3(3)];
    Q1Y1 = [Q1x2(1),Q1x4(1)];
    Q1Y2 = [Q1x2(2),Q1x4(2)];
    Q1Y3 = [Q1x2(3),Q1x4(3)];
    
    
    eb1 = [ux1,vx1,wx1];
    eb2 = [uy1,vy1,wy1];
    eb3 = [uz1,vz1,wz1];
    temp = 0:(2*pi)/15 :2*pi;
    r0 = 0.04;
    x1r1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x1;
    x1r2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x2;
    x1r3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x3;
    x1r4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x4;
    
    % Create legend
    if(ShowLegend)
        legend('Origin', 'X', 'Y', 'Z');
    end
    lineQ1Xhandle = line('XData',Q1X1,'YData',Q1X2,'ZData',Q1X3,'LineWidth',2,'Color',[1 0 0]);
    lineQ1Yhandle = line('XData',Q1Y1,'YData',Q1Y2,'ZData',Q1Y3,'LineWidth',2,'Color',[1 0 0]);
    Q1circle1handle = line('XData',x1r1(:,1),...
        'YData',x1r1(:,2),...
        'ZData',x1r1(:,3),...
        'LineWidth',2,'Color',[1 0 0]);
    
    Q1circle2handle = line('XData',x1r2(:,1),...
    'YData',x1r2(:,2),...
    'ZData',x1r2(:,3),...
    'LineWidth',2,'Color',[1 0 0]);

    Q1circle3handle = line('XData',x1r3(:,1),...
    'YData',x1r3(:,2),...
    'ZData',x1r3(:,3),...
    'LineWidth',2,'Color',[1 0 0]);

    Q1circle4handle = line('XData',x1r4(:,1),...
    'YData',x1r4(:,2),...
    'ZData',x1r4(:,3),...
    'LineWidth',2,'Color',[1 0 0]);
    
    
    
    org2Handle = plot3(oxQ2, oyQ2, ozQ2, 'k.'); % mass center of the quadcopter1
    quivX2handle = quiver3(oxQ2, oyQ2, ozQ2, ux2, vx2, wx2,  'r', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    quivY2handle = quiver3(oxQ2, oyQ2, ozQ2, uy2, vy2, wy2,  'g', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    quivZ2handle = quiver3(oxQ2, oyQ2, ozQ2, uz2, vz2, wz2,  'b', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    len2 = 0.34; % distance from the center of rotor 1 to rotor 3; 
    Q2x1 = [oxQ2,oyQ2,ozQ2] + 0.5*len2*[ux2,vx2,wx2];
    Q2x3 = [oxQ2,oyQ2,ozQ2] - 0.5*len2*[ux2,vx2,wx2];
    Q2x2 = [oxQ2,oyQ2,ozQ2] + 0.5*len2*[uy2,vy2,wy2];
    Q2x4 = [oxQ2,oyQ2,ozQ2] - 0.5*len2*[uy2,vy2,wy2];
    Q2X1 = [Q2x1(1),Q2x3(1)];
    Q2X2 = [Q2x1(2),Q2x3(2)];
    Q2X3 = [Q2x1(3),Q2x3(3)];
    Q2Y1 = [Q2x2(1),Q2x4(1)];
    Q2Y2 = [Q2x2(2),Q2x4(2)];
    Q2Y3 = [Q2x2(3),Q2x4(3)];
    
    
    eb1 = [ux2,vx2,wx2];
    eb2 = [uy2,vy2,wy2];
    eb3 = [uz2,vz2,wz2];
    temp = 0:(2*pi)/15 :2*pi;
    r0 = 0.04;
    x2r1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x1;
    x2r2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x2;
    x2r3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x3;
    x2r4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x4;
    
    % Create legend
    if(ShowLegend)
        legend('Origin', 'X', 'Y', 'Z');
    end
    lineQ2Xhandle = line('XData',Q2X1,'YData',Q2X2,'ZData',Q2X3,'LineWidth',2,'Color',[1 0 0]);
    lineQ2Yhandle = line('XData',Q2Y1,'YData',Q2Y2,'ZData',Q2Y3,'LineWidth',2,'Color',[1 0 0]);
    Q2circle1handle = line('XData',x2r1(:,1),...
        'YData',x2r1(:,2),...
        'ZData',x2r1(:,3),...
        'LineWidth',2,'Color','m');
    
    Q2circle2handle = line('XData',x2r2(:,1),...
    'YData',x2r2(:,2),...
    'ZData',x2r2(:,3),...
    'LineWidth',2,'Color','m');

    Q2circle3handle = line('XData',x2r3(:,1),...
    'YData',x2r3(:,2),...
    'ZData',x2r3(:,3),...
    'LineWidth',2,'Color','m');

    Q2circle4handle = line('XData',x2r4(:,1),...
    'YData',x2r4(:,2),...
    'ZData',x2r4(:,3),...
    'LineWidth',2,'Color','m');
    
    PayOrgHandle = plot3(oxP,oyP,ozP,'g.');% mass center of the payload
    
%     line1XLhandle = line('XData',[2-0.45,2-0.45],'YData',[3-0.0,3-0.0],'ZData',0.8+[-0.2,1.3],'LineWidth',4,'Color',[128 42 42]/256.0);
    line1XRhandle = line('XData',[2+0.5,2+0.5],'YData',[3+0.0,3+0.0],'ZData',0.8+[-0.2,1.3],'LineWidth',4,'Color',[128 42 42]/256.0);
    
    [sx,sy,sz] = sphere;
    SPx = 0.01*sx + oxP;
    SPy = 0.01*sy + oyP;
    SPz = 0.01*sz + ozP;
    payloadhandle = surf(SPx,SPy,SPz,'FaceColor',[1 0 0],'EdgeColor','none');
    cable1handle = line('Xdata',[oxQ1,oxP],...
        'Ydata',[oyQ1,oyP],...
        'Zdata',[ozQ1,ozP],...
        'LineWidth',2,'Color','k');
    cable2handle = line('Xdata',[oxQ2,oxP],...
        'Ydata',[oyQ2,oyP],...
        'Zdata',[ozQ2,ozP],...
        'LineWidth',2,'Color','k');
    % Window Frame1 (x:2.5,y:0)
%     line1XLhandle = line('XData',[3-0.4,3-0.4],'YData',[0,0],'ZData',0+[-0.2,1.3],'LineWidth',4,'Color',[128 42 42]/256.0);
%     line1XRhandle = line('XData',[3+0.4,3+0.4],'YData',[0,0],'ZData',0+[-0.2,1.3],'LineWidth',4,'Color',[128 42 42]/256.0);
%     line1ZUhandle = line('XData',[3-0.4,3+0.4],'YData',[0,0],'ZData',0+[1.3,1.3],'LineWidth',4,'Color',[128 42 42]/256.0);
%     line1ZDhandle = line('XData',[3-0.4,3+0.4],'YData',[0,0],'ZData',0+[-0.2,-0.2],'LineWidth',4,'Color',[128 42 42]/256.0);
    
    line1ZUhandle = line('XData',[2-0.5,2+0.5],'YData',[3-0.0,3+0.0],'ZData',0.8+[1.3,1.3],'LineWidth',4,'Color',[128 42 42]/256.0);
    
    % Set initial limits
    Xlim = [-2,4.0];
    Ylim = [0,7.5];
    Zlim = [-0,2.5];
    
%     Xlim = [xP(1)-AxisLength xP(1)+AxisLength] * LimitRatio;
%     Ylim = [yP(1)-AxisLength yP(1)+AxisLength] * LimitRatio;
%     Zlim = [zP(1)-AxisLength zP(1)+AxisLength] * LimitRatio;
    set(gca, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim);
    
        % Set initial view
    view(View(1, :));
%%
    for i = 1:numPlotSamples
        % Update graph title
        if(strcmp(Title, ''))
            titleText = sprintf('Sample %i of %i', 1+((i-1)*SamplePlotFreq), numSamples);
        else
            titleText = strcat(Title, ' (', sprintf('Sample %i of %i', 1+((i-1)*SamplePlotFreq), numSamples), ')');
        end
        title(titleText);
        
        % Plot body x y z axes
        if(strcmp(Trail, 'DotsOnly') || strcmp(Trail, 'All'))
            xQ1(1:i) = pQ1(1:i,1);
            yQ1(1:i) = pQ1(1:i,2);
            zQ1(1:i) = pQ1(1:i,3);
            oxQ1 = pQ1(i,1);% current mass center of Quadcopter 1
            oyQ1 = pQ1(i,2);
            ozQ1 = pQ1(i,3);
            
            xQ2(1:i) = pQ2(1:i,1);
            yQ2(1:i) = pQ2(1:i,2);
            zQ2(1:i) = pQ2(1:i,3);
            oxQ2 = pQ2(i,1);% current mass center of Quadcopter 2
            oyQ2 = pQ2(i,2);
            ozQ2 = pQ2(i,3);
            
            xP(1:i) = pP(1:i,1);
            yP(1:i) = pP(1:i,2);
            zP(1:i) = pP(1:i,3);
            oxP = pP(i,1);
            oyP = pP(i,2);
            ozP = pP(i,3);
            
            SPx = 0.05*sx + oxP;
            SPy = 0.05*sy + oyP;
            SPz = 0.05*sz + ozP;
        else
            xQ1 = pQ1(i,1);
            yQ1 = pQ1(i,2);
            zQ1 = pQ1(i,3);
            oxQ1 = pQ1(i,1);
            oyQ1 = pQ1(i,2);
            ozQ1 = pQ1(i,3);
            
            xQ2 = pQ2(i,1);
            yQ2 = pQ2(i,2);
            zQ2 = pQ2(i,3);
            oxQ2 = pQ2(i,1);
            oyQ2 = pQ2(i,2);
            ozQ2 = pQ2(i,3);
            
            xP = pP(i,1);
            yP = pP(i,2);
            zP = pP(i,3);
            oxP = pP(i,1);
            oyP = pP(i,2);
            ozP = pP(i,3);
            
            SPx = 0.01*sx + oxP;
            SPy = 0.01*sy + oyP;
            SPz = 0.01*sz + ozP;
            
        end
        if(strcmp(Trail, 'All'))
            ox1(1:i) = pQ1(1:i,1);
            oy1(1:i) = pQ1(1:i,2);
            oz1(1:i) = pQ1(1:i,3);
            ux1(1:i) = R1(1,1,1:i);
            vx1(1:i) = R1(2,1,1:i);
            wx1(1:i) = R1(3,1,1:i);
            uy1(1:i) = R1(1,2,1:i);
            vy1(1:i) = R1(2,2,1:i);
            wy1(1:i) = R1(3,2,1:i);
            uz1(1:i) = R1(1,3,1:i);
            vz1(1:i) = R1(2,3,1:i);
            wz1(1:i) = R1(3,3,1:i);
            
            ox2(1:i) = pQ2(1:i,1);
            oy2(1:i) = pQ2(1:i,2);
            oz2(1:i) = pQ2(1:i,3);
            ux2(1:i) = R2(1,1,1:i);
            vx2(1:i) = R2(2,1,1:i);
            wx2(1:i) = R2(3,1,1:i);
            uy2(1:i) = R2(1,2,1:i);
            vy2(1:i) = R2(2,2,1:i);
            wy2(1:i) = R2(3,2,1:i);
            uz2(1:i) = R2(1,3,1:i);
            vz2(1:i) = R2(2,3,1:i);
            wz2(1:i) = R2(3,3,1:i);
        else
            ox1 = pQ1(i,1);
            oy1 = pQ1(i,2);
            oz1 = pQ1(i,3);
            ux1 = R1(1,1,i);
            vx1 = R1(2,1,i);
            wx1 = R1(3,1,i);
            uy1 = R1(1,2,i);
            vy1 = R1(2,2,i);
            wy1 = R1(3,2,i);
            uz1 = R1(1,3,i);
            vz1 = R1(2,3,i);
            wz1 = R1(3,3,i);
            
            ox2 = pQ2(i,1);
            oy2 = pQ2(i,2);
            oz2 = pQ2(i,3);
            ux2 = R2(1,1,i);
            vx2 = R2(2,1,i);
            wx2 = R2(3,1,i);
            uy2 = R2(1,2,i);
            vy2 = R2(2,2,i);
            wy2 = R2(3,2,i);
            uz2 = R2(1,3,i);
            vz2 = R2(2,3,i);
            wz2 = R2(3,3,i);
        end
        % Update for Q1
        Q1x1 = [oxQ1,oyQ1,ozQ1] + 0.5*len1*[ux1,vx1,wx1];
        Q1x3 = [oxQ1,oyQ1,ozQ1] - 0.5*len1*[ux1,vx1,wx1];
        Q1x2 = [oxQ1,oyQ1,ozQ1] + 0.5*len1*[uy1,vy1,wy1];
        Q1x4 = [oxQ1,oyQ1,ozQ1] - 0.5*len1*[uy1,vy1,wy1];

        Q1X1 = [Q1x1(1),Q1x3(1)];
        Q1X2 = [Q1x1(2),Q1x3(2)];
        Q1X3 = [Q1x1(3),Q1x3(3)];
        Q1Y1 = [Q1x2(1),Q1x4(1)];
        Q1Y2 = [Q1x2(2),Q1x4(2)];
        Q1Y3 = [Q1x2(3),Q1x4(3)];
        
        eb1 = [ux1,vx1,wx1];
        eb2 = [uy1,vy1,wy1];
        eb3 = [uz1,vz1,wz1];
        temp = 0:(2*pi)/15 :2*pi;
        r0 = 0.04;
        x1r1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x1;
        x1r2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x2;
        x1r3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x3;
        x1r4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q1x4;
        
        % Update for Q2      
        Q2x1 = [oxQ2,oyQ2,ozQ2] + 0.5*len2*[ux2,vx2,wx2];
        Q2x3 = [oxQ2,oyQ2,ozQ2] - 0.5*len2*[ux2,vx2,wx2];
        Q2x2 = [oxQ2,oyQ2,ozQ2] + 0.5*len2*[uy2,vy2,wy2];
        Q2x4 = [oxQ2,oyQ2,ozQ2] - 0.5*len2*[uy2,vy2,wy2];

        Q2X1 = [Q2x1(1),Q2x3(1)];
        Q2X2 = [Q2x1(2),Q2x3(2)];
        Q2X3 = [Q2x1(3),Q2x3(3)];
        Q2Y1 = [Q2x2(1),Q2x4(1)];
        Q2Y2 = [Q2x2(2),Q2x4(2)];
        Q2Y3 = [Q2x2(3),Q2x4(3)];
        
        eb1 = [ux2,vx2,wx2];
        eb2 = [uy2,vy2,wy2];
        eb3 = [uz2,vz2,wz2];
        temp = 0:(2*pi)/15 :2*pi;
        r0 = 0.04;
        x2r1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x1;
        x2r2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x2;
        x2r3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x3;
        x2r4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*Q2x4;
        %
        
        ratio = 0.4;
%         set(org1Handle, 'xdata', xQ1, 'ydata', yQ1, 'zdata', zQ1,'Color','b','MarkerSize',3);% Q1 mass center
        set(quivX1handle, 'xdata', oxQ1, 'ydata', oyQ1, 'zdata', ozQ1,'udata', ratio*ux1, 'vdata', ratio*vx1, 'wdata', ratio*wx1);
        set(quivY1handle, 'xdata', oxQ1, 'ydata', oyQ1, 'zdata', ozQ1,'udata', ratio*uy1, 'vdata', ratio*vy1, 'wdata', ratio*wy1);
        set(quivZ1handle, 'xdata', oxQ1, 'ydata', oyQ1, 'zdata', ozQ1,'udata', ratio*uz1, 'vdata', ratio*vz1, 'wdata', ratio*wz1);
        set(lineQ1Xhandle,'XData',Q1X1,'YData',Q1X2,'ZData',Q1X3);
        set(lineQ1Yhandle,'XData',Q1Y1,'YData',Q1Y2,'ZData',Q1Y3);
        set(Q1circle1handle,'XData',x1r1(:,1),'YData',x1r1(:,2),'ZData',x1r1(:,3));
        set(Q1circle2handle,'XData',x1r2(:,1),'YData',x1r2(:,2),'ZData',x1r2(:,3));
        set(Q1circle3handle,'XData',x1r3(:,1),'YData',x1r3(:,2),'ZData',x1r3(:,3));
        set(Q1circle4handle,'XData',x1r4(:,1),'YData',x1r4(:,2),'ZData',x1r4(:,3));
        
%         set(org2Handle, 'xdata', xQ2, 'ydata', yQ2, 'zdata', zQ2,'Color','c','MarkerSize',3);% Q2 mass center
        set(quivX2handle, 'xdata', oxQ2, 'ydata', oyQ2, 'zdata', ozQ2,'udata', ratio*ux2, 'vdata', ratio*vx2, 'wdata', ratio*wx2);
        set(quivY2handle, 'xdata', oxQ2, 'ydata', oyQ2, 'zdata', ozQ2,'udata', ratio*uy2, 'vdata', ratio*vy2, 'wdata', ratio*wy2);
        set(quivZ2handle, 'xdata', oxQ2, 'ydata', oyQ2, 'zdata', ozQ2,'udata', ratio*uz2, 'vdata', ratio*vz2, 'wdata', ratio*wz2);
        set(lineQ2Xhandle,'XData',Q2X1,'YData',Q2X2,'ZData',Q2X3);
        set(lineQ2Yhandle,'XData',Q2Y1,'YData',Q2Y2,'ZData',Q2Y3);
        set(Q2circle1handle,'XData',x2r1(:,1),'YData',x2r1(:,2),'ZData',x2r1(:,3));
        set(Q2circle2handle,'XData',x2r2(:,1),'YData',x2r2(:,2),'ZData',x2r2(:,3));
        set(Q2circle3handle,'XData',x2r3(:,1),'YData',x2r3(:,2),'ZData',x2r3(:,3));
        set(Q2circle4handle,'XData',x2r4(:,1),'YData',x2r4(:,2),'ZData',x2r4(:,3));
        
        set(PayOrgHandle,'xdata', xP, 'ydata', yP, 'zdata', zP,'Color','g','LineWidth',0.1);
        set(payloadhandle,'XData',SPx,'YData',SPy,'ZData',SPz,'FaceColor',[1 0 0],'EdgeColor','none');    
        set(cable1handle, 'XData',[oxQ1,oxP],'YData',[oyQ1,oyP],'ZData',[ozQ1,ozP]);
        set(cable2handle, 'XData',[oxQ2,oxP],'YData',[oyQ2,oyP],'ZData',[ozQ2,ozP]);
        
        if(isFixView)
             Xlim(1) = pP(i,1) - LimitRatio*AxisLength;
             Ylim(1) = pP(i,2) - LimitRatio*AxisLength;
             Zlim(1) = pP(i,3) - LimitRatio*AxisLength;
             Xlim(2) = pP(i,1) + LimitRatio*AxisLength;
             Ylim(2) = pP(i,2) + LimitRatio*AxisLength;
             Zlim(2) = pP(i,3) + LimitRatio*AxisLength;
             set(gca, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim);
        else
    %       % Adjust axes for snug fit and draw
            axisLimChanged = false;
            if((pP(i,1) - AxisLength) < Xlim(1)), Xlim(1) = pP(i,1) - LimitRatio*AxisLength; axisLimChanged = true; end
            if((pP(i,2) - AxisLength) < Ylim(1)), Ylim(1) = pP(i,2) - LimitRatio*AxisLength; axisLimChanged = true; end
            if((pP(i,3) - AxisLength) < Zlim(1)), Zlim(1) = pP(i,3) - LimitRatio*AxisLength; axisLimChanged = true; end
            if((pP(i,1) + AxisLength) > Xlim(2)), Xlim(2) = pP(i,1) + LimitRatio*AxisLength; axisLimChanged = true; end
            if((pP(i,2) + AxisLength) > Ylim(2)), Ylim(2) = pP(i,2) + LimitRatio*AxisLength; axisLimChanged = true; end
            if((pP(i,3) + AxisLength) > Zlim(2)), Zlim(2) = pP(i,3) + LimitRatio*AxisLength; axisLimChanged = true; end
            if(axisLimChanged), set(gca, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim); end
        end
        drawnow;
% 
        %%% Adjust view
        if(CreateAVI)
%             if(numel(View) > 2)
%                 view(View(i, :));
%             end
        end
% 
        % Add frame to AVI object
        if(~isempty(aviobj))
            frame = getframe(fig);
            writeVideo(aviobj,frame);
%             aviobj = addframe(aviobj, frame);
        end

    end

end

