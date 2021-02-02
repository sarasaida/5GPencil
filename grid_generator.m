%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% GRID AND LAYOUT GENERATOR %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%To define the cells layout and to create the meshgrid for the generation
%of the 'Deployment spot'(candidate points where the beam will be directed)
%and 'Measure spot' (to measure the power density and the electric field)

% Authors: Simone Rossetti, Sara Saida
% Date: 21.01.2021
% Version: 1.1
%% Parameters definition

close all
clear all
clf
clear

tic

r=100; %cell radius
spacing=1; %spacing of the meshgrid in [m]
num_user=64; %number of user per sector (max 64)
measure=0; %1 for measure spot, 0 for deploy spot 
plot_on=0; %1 to enable plot, 0 to disable plot

%% Layout Definition

%CREATION OF THE HEXAGONAL GRID centered in (0,0) 
%location of the centers of the 6 neighbouring hex BSs
BS=[0 0]; %center of serving BS

for theta=(30):(60):(330)
    x= r*sqrt(3)* cosd(theta);
    y= r*sqrt(3)* sind(theta);
    BS = [BS; x y];
end

%creation of the serving hex BS and neighbouring hex BSs
for i=1:length(BS)
    xc = BS(i,1);
    yc = BS(i,2);
    xn = [(r+xc) (r/2+xc) (-r/2+xc) (-r+xc) (-r/2+xc) (r/2+xc) (r+xc)];
    yn = [yc (r*sqrt(3)/2+yc) (r*sqrt(3)/2+yc) yc (-r*sqrt(3)/2+yc) (-r*sqrt(3)/2+yc) yc];
    cell{i}=polyshape(xn,yn);
end

for i=1:length(BS)
    xc = BS(i,1);
    yc = BS(i,2);
   
    %sector1
    x1 = [xc+3*r/4 xc+r/2 xc-r/2 xc-3*r/4 xc];
    y1 = [yc+r*sqrt(3)/4 yc+r*sqrt(3)/2 yc+r*sqrt(3)/2 yc+r*sqrt(3)/4 yc];
    %sector2
    x2 =[xc-3*r/4 xc-r xc-r/2 xc xc];
    y2 =[yc+r*sqrt(3)/4 yc yc-r*sqrt(3)/2 yc-r*sqrt(3)/2 yc];
    %sector3
    x3 =[xc xc+r/2 xc+r xc+3*r/4 xc];
    y3 =[yc-r*sqrt(3)/2 yc-r*sqrt(3)/2 yc yc+r*sqrt(3)/4 yc];
    
    %area not accessible to the user, to be subtracted from the sector
    radius = 10;  
    th = 0:pi/50:2*pi;
    xunit = radius * cos(th) + xc;
    yunit = radius * sin(th) + yc;
    poly0 = polyshape(xunit(1:end-1),yunit(1:end-1));
    
    %create full sectors
    poly1 = polyshape(x1,y1);
    poly2 = polyshape(x2,y2);
    poly3 = polyshape(x3,y3);
    
    %usable area of the sectors
    sector{i,1} = subtract(poly1, poly0);
    sector{i,2} = subtract(poly2, poly0);
    sector{i,3} = subtract(poly3, poly0);
    
    
    %sector 1 right side
    x1_rx = [xc+3*r/4 xc+r/2 xc xc ];
    y1_rx = [yc+r*sqrt(3)/4 yc+r*sqrt(3)/2 yc+r*sqrt(3)/2 yc];
    %sector 1 left side
    x1_lx =[xc  xc-r/2 xc-3*r/4 xc ];
    y1_lx =[yc+r*sqrt(3)/2 yc+r*sqrt(3)/2 yc+r*sqrt(3)/4 yc];
    %create right side of the sectors
    secside_r{i,1} = polyshape(x1_rx,y1_rx);
    %create left side of the sectors
    secside_l{i,1} = polyshape(x1_lx,y1_lx);
    
    %sector 2 right side
    x2_rx = [xc-3*r/4 xc-r xc-3*r/4 xc ];
    y2_rx = [yc+r*sqrt(3)/4 yc yc-r*sqrt(3)/4 yc ];
    %sector 2 left side
    x2_lx =[xc-3*r/4  xc-r/2 xc xc ];
    y2_lx =[yc-r*sqrt(3)/4 yc-r*sqrt(3)/2 yc-r*sqrt(3)/2 yc];
    %create right side of the sector
    secside_r{i,2} = polyshape(x2_rx,y2_rx);
    %create left side of the sector
    secside_l{i,2} = polyshape(x2_lx,y2_lx);
    
    %sector 3 right side
    x3_rx = [xc xc+r/2 xc+3*r/4 xc ];
    y3_rx = [yc-r*sqrt(3)/2 yc-r*sqrt(3)/2 yc-r*sqrt(3)/4 yc ];
    %sector 3 left side
    x3_lx =[xc+3*r/4  xc+r xc+3*r/4 xc ];
    y3_lx =[yc-r*sqrt(3)/4 yc yc+r*sqrt(3)/4 yc];
    %create right side of the sectors
    secside_r{i,3} = polyshape(x3_rx,y3_rx);
    %create left side of the sectors
    secside_l{i,3} = polyshape(x3_lx,y3_lx);
    
end

if plot_on==1
    for j=1:length(cell)
        for i=1:size(sector,2)
            plot(sector{j,i});
            axis equal;
            hold on;
        end
    end
end
save('layout.mat','BS', 'cell', 'sector', 'secside_r', 'secside_l');

%% Meshgrid Generation

if measure==1
    %define the ends of the grid
    min_x = -r;
    max_x = r;

    min_y = -sqrt(3)*r/2;
    max_y = sqrt(3)*r/2;

    %create meshgrid 
    x = min_x:spacing:max_x;
    y = min_y:spacing:max_y;
    [X,Y] = meshgrid(x,y);
    total_spots = [X(:), Y(:)];

    %% Check area of ​​belonging of the useful spots

    spots_cat{length(cell),size(sector,2)}=[]; %spots cataloged by area of ​​belonging (cellxsector)

    for i=1:length(total_spots)
        spot=[total_spots(i,1),total_spots(i,2)];
        for j=1:length(cell)
            flag = isinterior(cell{j},spot); %check if it is a useful spot
            if flag == 1 %spot inside the cell
                flag1 = isinterior(sector{j,1},spot); %check sector 1
                flag2 = isinterior(sector{j,2},spot); %check sector 2
                flag3 = isinterior(sector{j,3},spot); %check sector 3

                if flag1 == 1 %spot inside sector 1
                    flagr = isinterior(secside_r{j,1},spot); 
                    if flagr ==1 %spot in right side
                        spots_cat{j,1} = [spots_cat{j,1}; spot,1]; %coord of spot, 1 if right 0 if left
    %                     useful_spots{end+1}=[spot, j, 1, 1]; %[x, y, num_cell, num_sec, rx (1) or lx(0)]
                    else %spot in left side
                        spots_cat{j,1} = [spots_cat{j,1}; spot,0];
                    end


                elseif flag2 == 1 %spot inside sector 2
                    flagr = isinterior(secside_r{j,2},spot); 
                    if flagr ==1 %spot in right side
                        spots_cat{j,2} = [spots_cat{j,2}; spot, 1];
                    else %spot in left side
                        spots_cat{j,2} = [spots_cat{j,2}; spot, 0];
                    end

                elseif flag3 == 1 %spot inside sector 3
                    flagr = isinterior(secside_r{j,3},spot); 
                    if flagr ==1 %spot in right side
                        spots_cat{j,3} = [spots_cat{j,3}; spot, 1];
                    else %spot in left side
                        spots_cat{j,3} = [spots_cat{j,3}; spot, 0];
                    end
                end
            end
        end
    end

    save('useful_spots.mat','spots_cat');

    %% Measure Spot selection
    %Select a subset of points or the complete set of points inside the cell of serving
    %BS to measure the electric field

    measure_spot{1,size(spots_cat,2)}=[];

    for i=1:size(spots_cat,2)
        for j=1:length(spots_cat{1,i})
            measure_spot{1,i} = [measure_spot{1,i}; spots_cat{1,i}(j,:)];
        end
    end

    save('measure.mat','measure_spot');
    if plot_on==1
        for i=1:length(measure_spot)
            for k=1:length(measure_spot{1,i})
                plot(measure_spot{1,i}(k,1),measure_spot{1,i}(k,2),'k.');
                hold on
            end
        end
    end
    
elseif measure==0
    %% Deploy Spot selection
    % Pseudorandomly generate deploy spot exploiting polar coordinates 
    
    deploy_spot{length(cell),size(sector,2)}=[];

    for j=1:length(cell)
        for i=1:size(sector,2)
            for k=1:num_user
                xc = BS(j,1);
                yc = BS(j,2);
                c=0;
                if i==1 %sector 1
                    theta=30:0.1:150;
                elseif i==2 %sector 2
                    theta=150:0.1:270;
                elseif i==3 %sector 3
                    theta=-90:0.1:30;
                end
                while c<1
                    c=1;
                    key = randi([1, length(theta)]);
                    r1=((r-radius)*rand(1,1))+radius; %select r1 between 0-90 and add 10
                    theta1=theta(key);
                    x1=xc+r1*cosd(theta1); % convert to cartesian
                    y1=yc+r1*sind(theta1);
                    spot=[x1,y1];
                    if isinterior(secside_r{j,i},spot)==1 %sector side rx
                        deploy_spot{j,i} = [deploy_spot{j,i}; spot,1];
                    elseif isinterior(secside_l{j,i},spot)==1 %sector side lx
                        deploy_spot{j,i} = [deploy_spot{j,i}; spot,0];
                    else
                        c=0; %if falls outside sector then regenerate
                    end
                end
            end
        end
    end
    
    if plot_on==1
        for j=1:length(cell)
            for i=1:size(deploy_spot,2)
                for k=1:length(deploy_spot{j,i})
                    plot(deploy_spot{j,i}(k,1),deploy_spot{j,i}(k,2),'b*');
                    axis equal;
                    hold on;

                end
            end
        end
    end
    save('deploy.mat','deploy_spot');    
    
    
end
    
toc
    
    



