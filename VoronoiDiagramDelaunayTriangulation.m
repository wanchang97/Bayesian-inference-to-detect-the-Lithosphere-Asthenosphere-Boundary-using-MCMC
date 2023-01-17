close all; clear all; clc
%% Voronoi diagram and delaunayTriangulation
lx = 4;
ly = 4;
Cx_true = lx*[1/5 1/3 2/5 7/8 1/7]';
Cy_true = ly*[3/4 1/2 1/3 3/4 2/5]';
%C_true = [1/5*lx 1/4*ly;
%          1/3*lx 1/2*ly;
%          2/5*lx 1/3*ly;
%          7/8*lx 3/4*ly;
%          1/7*lx 1/5*ly];
%DT = delaunay(C_true);
%plot(C_true(:,1),C_true(:,2),'*r')
DT = delaunayTriangulation(Cx_true,Cy_true);
IC = incenter(DT);

figure
subplot(2,2,1);
voronoi(Cx_true,Cy_true)
hold on
plot(Cx_true,Cy_true,'*r')
xlim([0,lx]);
ylim([0,ly]);
xlabel('X');
ylabel('Y');
title('Voronoi diagram')
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
axis equal

subplot(2,2,2);
triplot(DT,'g-.','LineWidth',2)
hold on
plot(Cx_true,Cy_true,'*r')
xlim([0,lx]);
ylim([0,ly]);
xlabel('X');
ylabel('Y');
title('Delaunay triangulation')
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
axis equal

subplot(2,2,[3,4]);
triplot(DT,'g-.','LineWidth',2)
hold on
plot(Cx_true,Cy_true,'*r')
voronoi(Cx_true,Cy_true)
plot(Cx_true,Cy_true,'*r')
xlim([0,lx]);
ylim([0,ly]);
xlabel('X');
ylabel('Y');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
axis equal
title('Both Delaunay triangulation and Voronoi diagram')
saveas(gcf,'VoronoiDiagramDelaunayTriangulation.png')