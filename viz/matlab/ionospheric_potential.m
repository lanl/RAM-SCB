% Plots ionospheric and equatorial (mapped along SC B-field lines) electric potential 
% Optionally, interpolates equatorial potential onto RAM grid and
% writes file

% S. Zaharia, 2010

% function ionospheric_potential(hour)

close all;
clear all;

hour = 12.08; % Hour of simulation where results (SCB) are plotted
raminterp = true; % If true, interpolates/writes values on RAM grid
                  % for all times in the file
plotram = false; % If true, plots values on RAM grid for all times
                 % in the file (use in moderation; can by large # of figures)

% RAM grid
nr = 20;
nphi = 25;
radOut = 6.75;

strHour = num2str(hour, '%05.2f')

HOME = getenv('HOME');

% Location of ionospheric.nc file
PREFIX = strcat(HOME,'/3Dcode/RAM-SCB_SAND/RAM_SCB/run_scb_geo/IM/output_scb/');
filename = strcat(PREFIX, 'ionospheric_potential.nc')

% Open NETCDF file and read variables
ncid = netcdf.open(filename, 'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'time');
time = netcdf.getVar(ncid,varid,'double')
% Find time index iTime closest to hour
iTime = find(time>=hour,1)
varid = netcdf.inqVarID(ncid,'xIono');
xIonoFull = netcdf.getVar(ncid,varid,'double');
xIono = squeeze(xIonoFull(:,:,iTime));
varid = netcdf.inqVarID(ncid,'yIono')
yIonoFull = netcdf.getVar(ncid,varid,'double');
yIono = squeeze(yIonoFull(:,:,iTime));
varid = netcdf.inqVarID(ncid,'xEq')
xEqt = netcdf.getVar(ncid,varid,'double');
xEq = squeeze(xEqt(:,:,iTime));
varid = netcdf.inqVarID(ncid,'yEq')
yEqt = netcdf.getVar(ncid,varid,'double');
yEq = squeeze(yEqt(:,:,iTime));
varid = netcdf.inqVarID(ncid,'PhiIono')
PhiIonot = netcdf.getVar(ncid,varid,'double');
PhiIono = squeeze(PhiIonot(:,:,iTime));
% Close NETCDF file
netcdf.close(ncid);


PhiIono = 1e-3*PhiIono; % in kV

nalpha = size(xIono,1)
nbeta = size(xIono,2)
halfbeta = fix(nbeta/2) + 1;

phimin = min(min(PhiIono))
phiMax = max(max(PhiIono))

max(max(PhiIono))-min(min(PhiIono))

figure(1);
pcolor(xEq,yEq,PhiIono);
hold on;

for i = phimin:5:phiMax
  v = [i, i];
  c = contour(xEq, yEq, PhiIono, v, 'LineWidth', 2, 'Color', 'black');
  clabel(c);
end

h = colorbar;
set(h,'FontName','Helvetica','FontSize',16, 'FontWeight', ...
      'bold');
 set(h, 'YLimMode', 'auto', 'YTickMode', 'auto');
 set(h, 'TickDir', 'out');
% Set axes equal, title and reorient the figure
axis equal;
set(gca,'FontName','Helvetica','FontSize',16, 'FontWeight','bold');
title(strcat('Convective potential mapped onto eq. plane at T = ',strHour,'(kV)'));
shading flat;
xlabel('X');
ylabel('Y');
view(180,90);

t = 0:0.01:2*pi;
x1 = cos(t);
y1 = sin(t);
plot(x1,y1, 'k');
t2 = pi/2:0.01:3*pi/2;
for i = 1:100
  x2 = i/100*cos(t2);
  y2 = i/100*sin(t2);
  z = ones(size(x2,2),size(y2,2));
  plot(x2, y2, 'k');
end

figure(2);
pcolor(xEq, yEq, PhiIono);
% set(gca, 'clim', [-30 0]);
h = colorbar;
set(h,'FontName','Helvetica','FontSize',16, 'FontWeight', ...
      'bold');
set(h, 'YLimMode', 'auto', 'YTickMode', 'auto');
set(h, 'TickDir', 'out');
% Set axes equal, title and reorient the figure
axis equal;
set(gca,'FontName','Helvetica','FontSize',16, 'FontWeight','bold');
title(strcat('Convective potential in eq. plane at T = ',strHour,'(kV)'));
shading flat;
xlabel('X');
ylabel('Y');
view(180,90);
shading faceted;

hold on;
t = 0:0.01:2*pi;
x1 = cos(t);
y1 = sin(t);
plot(x1,y1, 'k');
t2 = pi/2:0.01:3*pi/2;
for i = 1:100
  x2 = i/100*cos(t2);
  y2 = i/100*sin(t2);
  z = ones(size(x2,2),size(y2,2));
  plot(x2, y2, 'k');
end

figure(3);
  pcolor(xIono, yIono, PhiIono);
  % set(gca, 'clim', [-30 30]);
  h = colorbar;
  set(h,'FontName','Helvetica','FontSize',16, 'FontWeight', ...
	  'bold');
  set(h, 'YLimMode', 'auto', 'YTickMode', 'auto');
  set(h, 'TickDir', 'out');
  % Set axes equal, title and reorient the figure
  axis equal;
  set(gca,'FontName','Helvetica','FontSize',16, 'FontWeight','bold');
  title(strcat('Convective potential on ionosphere at T = ',strHour,' (kV)'));
  shading flat;
  xlabel('X');
  ylabel('Y');
  view(-90,90);
  axis off;
  
  hold on;
  
  theta = 0:.01:3*pi;
  x = sin(10.*pi/180.) * sin(theta);
  y = sin(10.*pi/180.) * cos(theta);
  
  plot(x, y, 'k');
  % text(0.75,0.75, '80^\circ', 'fontsize', 13, 'FontWeight', 'normal');
  hold on;
  fac2 = sind(20.)/sind(10.)
  plot(fac2*x, fac2*y, 'k');
  % text(1.4, 1.4, '70^\circ', 'fontsize', 13, 'FontWeight', 'normal');
  
  fac3 = sind(30.)/sind(10.)
  plot(fac3*x, fac3*y, 'k');
  % text(2.15, 2.15, '60^\circ', 'fontsize', 13, 'FontWeight', 'normal');
  
  fac4 = sind(40.)/sind(10.)
  plot(fac4*x, fac4*y, 'k');
  % text(2.95, 2.95, '50^\circ', 'fontsize', 13, 'FontWeight', 'normal');

  set(gcf, 'Position', [172 144 1264 876]);

 if (raminterp)
   
  radRAM = zeros(nr,1);
  azimRAM = zeros(nphi,1);
  xRAM = zeros(nr,nphi);
  yRAM = zeros(nr,nphi);
  
  xscatter = zeros(nalpha*(nbeta-1),1);
  yscatter = zeros(nalpha*(nbeta-1),1);
  phiscatter = zeros(nalpha*(nbeta-1),1);
  
  for j = 1:nr
     radRAM(j) = 2.0 + (radOut-2.0) * (j-1.)/(nr-1.);
%    ! print*, 'hRAM: j1, radRAM(j1): ', j1, radRAM(j1)
  end
  
  for k = 1:nphi
     azimRAM(k) = 24. * (k-1.)/(nphi-1.);
  end
  
  for j = 1:nr
    for k = 1:nphi
      xRAM(j,k) = radRAM(j) .* cos(azimRAM(k)*2.*pi/24. - pi);
      yRAM(j,k) = radRAM(j) .* sin(azimRAM(k)*2.*pi/24. - pi);
    end
  end
 
  
  ntimes = size(time,1);
  
  for itime = 1:ntimes
  
    inum = 0;
    
    for j = 1:nalpha
      for k = 2:nbeta
        inum = inum+1;
        xscatter(inum) = xEqt(j,k,itime);
        yscatter(inum) = yEqt(j,k,itime);
        phiscatter(inum) = 1e-3*PhiIonot(j,k,itime);
      end
    end
    
    F = TriScatteredInterp(xscatter,yscatter,phiscatter, 'natural');
    phiRAM = F(xRAM,yRAM);
    
    strHour = num2str(time(itime), '%05.2f')
    
    if (plotram | itime == iTime)
      figure(4+itime);
      pcolor(xRAM,yRAM,phiRAM);
      set(gca, 'clim', [phimin phiMax]);
      h = colorbar;
      set(h,'FontName','Helvetica','FontSize',16, 'FontWeight', ...
            'bold');
      set(h, 'YLimMode', 'auto', 'YTickMode', 'auto');
      set(h, 'TickDir', 'out');
      
      axis equal;
      set(gca,'FontName','Helvetica','FontSize',16, 'FontWeight','bold');
      title(strcat('Convective potential on RAM grid at T = ',strHour,' (kV)'));
      shading flat;
      xlabel('X');
      ylabel('Y');
      view(180,90);
      shading faceted;
      
      hold on;
      t = 0:0.01:2*pi;
      x1 = cos(t);
      y1 = sin(t);
      plot(x1,y1, 'k');
      t2 = pi/2:0.01:3*pi/2;
      for i = 1:100
        x2 = i/100*cos(t2);
        y2 = i/100*sin(t2);
        z = ones(size(x2,2),size(y2,2));
        plot(x2, y2, 'k');
      end
    end
      
  fid = fopen(strcat('EPOT_SC_',strHour,'.txt'), 'w');
  fprintf(fid, '%s\n', 'Hour of simulation');
  fprintf(fid, '%10.3f\n', time(itime));
  fprintf(fid, '%s\n', 'Ionospheric potential mapped along SCB lines');
  fprintf(fid, '%s\n',  'L       PHI            Epot(kV)');
 
  for i = 1:nr
    for j = 1:nphi
      fprintf(fid, '%4.2f %10.6f %13.4f\n', radRAM(i), azimRAM(j)*2*pi/24., phiRAM(i,j));
    end
  end
  fclose(fid);
  
  end
  
 end

 figure;
 subplot1(1,2);
 subplot1(1);
 
 pcolor(xEq,yEq,PhiIono);
hold on;

for i = phimin:5:phiMax
  v = [i, i];
  c = contour(xEq, yEq, PhiIono, v, 'LineWidth', 2, 'Color', 'black');
  clabel(c);
end

% h = colorbar;
% set(h,'FontName','Helvetica','FontSize',16, 'FontWeight', ...
%      'bold');
% set(h, 'YLimMode', 'auto', 'YTickMode', 'auto');
% set(h, 'TickDir', 'out');
% Set axes equal, title and reorient the figure
axis image;
set(gca,'FontName','Helvetica','FontSize',16, 'FontWeight','bold');
title(strcat('Convective potential mapped onto eq. plane at T = ',strHour,'(kV)'));
shading flat;
xlabel('X');
ylabel('Y');
view(180,90);

t = 0:0.01:2*pi;
x1 = cos(t);
y1 = sin(t);
plot(x1,y1, 'k');
t2 = pi/2:0.01:3*pi/2;
for i = 1:100
  x2 = i/100*cos(t2);
  y2 = i/100*sin(t2);
  z = ones(size(x2,2),size(y2,2));
  plot(x2, y2, 'k');
end
pos1 = get(gca, 'Position'); 

subplot1(2);
 pcolor(xRAM,yRAM,phiRAM);
      set(gca, 'clim', [phimin phiMax]);
      h = colorbar;
      set(h,'FontName','Helvetica','FontSize',16, 'FontWeight', ...
            'bold');
      set(h, 'YLimMode', 'auto', 'YTickMode', 'auto');
      set(h, 'TickDir', 'out');
      
      axis image;
      set(gca,'FontName','Helvetica','FontSize',16, 'FontWeight','bold');
      title(strcat('Convective potential on RAM grid at T = ',strHour,' (kV)'));
      shading flat;
      xlabel('X');
      % ylabel('Y');
      view(180,90);
      shading faceted;
      
      hold on;
      t = 0:0.01:2*pi;
      x1 = cos(t);
      y1 = sin(t);
      plot(x1,y1, 'k');
      t2 = pi/2:0.01:3*pi/2;
      for i = 1:100
        x2 = i/100*cos(t2);
        y2 = i/100*sin(t2);
        z = ones(size(x2,2),size(y2,2));
        plot(x2, y2, 'k');
      end
      set(gca, 'Position', pos1+[0.425 0 0 0]);
      
   