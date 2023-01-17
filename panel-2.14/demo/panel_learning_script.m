% Panel learning script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanel2 basic use as a subplot
% a) Create a grid of panels
% b) Plot into each sub-panel

% a) create a NxN grid in gcf (this will create a figure, if none is open).
% You can pass the figure handle to the constructor if you need to attach
% the panel to a particular figure as p = panel(h_figure)
N = 2;
use_panel = 1;
clf % deletes all children of the current figure that have visible handles.
% Prepare
if use_panel
    p = panel();
    p.pack(N,N);
end
% b)
% plot into each panel in turn
for m = 1:N
    for n = 1: N
        % select on of the NxN grid of sub-panels
        if use_panel
            p(m,n).select();
        else
            subplot(N,N,m+(n-1)*N);
        end
        % plot some data
        plot(randn(100,1));
        % you can use all the usual calls
        xlabel('sample number');
        ylabel('data');
        axis([0 100 -3 3]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demopanel3
clf
% Pnel nesting
% a) Create a grid of panels
% b) Plot into three of the sub-panels.
% c) Create another grid in the fourth.
% d) Plot into each of these.
% a) Create a panel in gcf
% "p" is called the "root panel", which is the special panel whose parent
% is the figure window(usually), rather than another panel.
p = panel();
% pack a 2x2 grid of panels into it
p.pack(2,2);
% b) Plot into the first three panels
for m = 1:2
    for n = 1:2
        % skip the 2,2 panel
        if m == 2 && n == 2
            break
        end
        % select the panel (create an axis, and make that axis current)
        p(m,n).select();
        % plot some stuff
        plot(randn(100,1));
        xlabel('sample number');
        ylabel('data');
        axis([0 100 -3 3]);

    end
end
% c) Pack a further grid into p(2,2)
% all panels start as "uncommitted panels" (even the root panel). the first
% time we "select()" one, we commit it as an "axis panel". the first time
% we "pack()" one, we commit it as a "parent panel". once committed, it
% can't be changed into the other sort.
p(2,2).pack(2,3);
%% d) Plot into the six new sub-sub-panels
for m  =1:2
    for n = 1:3
        % select the panel - this commits it as an axis panel
        p(2,2,m,n).select();
        % plot some stuff
        plot(randn(100,1)); % normally distributed random number
        xlabel('sample number');
        ylabel('data');
        axis([0 100 -3 3]);
    end
end
%%
% note this alternative, equivalent, way to reference a sub-panel
p_22 = p(2,2);
% plot another bit of data into the six sub-sub-panels
for m = 1:2
    for n = 1:3
        % select the panel
        p_22(m,n).select();
        % plot more stuff
        hold on
        plot(randn(100,1)*0.3,'r');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanel4
clf
% panels can be any size
% a) Create an asymmetrical grid of panels
% b) Create another
% c) Use select('all') to load them all with axes
% d) Get handles to all the axes and modify them.
%% a)
% create a 2x2 grid in gcf with different fractionally-sized rows and
% columns. a row or column sized as "[]" with stretch to fill the remaining
% unassigned space
p = panel();
p.pack({1/3 []}, {1/3 []});
%% b) Pack a 2x3 grid into p(2,2). Note that we can pack by percentage as
% well as by fraction - the interpretation is just based on the size of the
% numbers we pass in(1 to 100 for percentage, or 0 to 1 for fraction).
p(2,2).pack({30 70},{20 20 []});
%% c) Use select('all') to quickly show the  layout you've achieved.This commits all uncommited panels as axis panels
% so they can't be parents anymore(i.e. they can't have more children
% pack()ed into them).
% this is no use once you've got organised- e.g. the first three demos
% don't use it. But it may help to see what you're doing in a starting
% stage
p.select('all');
%% d) Whilst we're here, we can get all the axes within a particular panel like this. There are three "groups" 
% associated with a panel: (fa)mily, (de)scendants, and (ch)ildre. see
% "help panel/descendants", for instance, to see who's in them. they're
% each useful in different circumstances. here we use (de)scendants.
h_axes = p.de.axis;
% so then we might want to set something on them.
set(h_axes,'color',[0 0 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanel5 Tools for finding your way around a layout
clf
% a) Recreate the complex layout from demopanel1
% b) Show threee tools that help to navigate a layout
figure(1)
clf
% create panel
p = panel();
% layout a variety of sub-panels
p.pack('h',{1/3 []})% 1x2 pack two columns. to pack columns instead
% of rows, we just pass "h" (horizontal) to pack() .
p(1).pack({2/3 []});% 2x1
p(1,1).pack(3,2);
p(2).pack(6,2);

% set margins
p.de.margin = 2;
p(1,1).marginbottom = 12;
p(2).marginleft = 20;
p.margin = [13 10 2 2];
% and some properties
p.fontsize = 8;
% if a layout gets complex, it can be tricky to find your way around it.
% it's quite natural once you get the hang, but there are three tools that
% will help you if you get lost: display(), identify() and show()
% identify() only works on axis panels. we haven't bothered plotting any
% dat, this time, so we'll use select('all') to commit all remaining
% uncommited panels as axis panels
p.select('all');
% display() the panel object at the prompt
% notice that most of the panels are called "Object"- this is because they
% are "object panels", which is the general name for axis panels (and
% that's because panels can contain other graphics objects as well as
% axes).
p

% use identify() every panel that is an axis panel has its axis wiped and
% replaced with the panel's reference. the one in the bottom right, for
% instance, is labelled "(2,6,2)", which means we can access it with
% p(2,6,2)
p.identify();

% use show() the selected panel is highlighted in red. shown works on
% parent panel as well 
p(2,6,2).show();
p(2).show();
% just to prove the point, let's now select one of the panels we've
% identified and plot something into it.
p(2,4,1).select();
plot(randn(100,1))
axis auto
%% b)
% date set 1
for m = 1:3
    for n = 1:2
        % prepare sample data
        t = (0:99) /100;
        s1 = sin(t*2*pi*m);
        s2 = sin(t*2*pi*n*2);
        % select axis - see data set 2 for an alternative way to access
        % sub-panels
        p(1,1,m,n).select();
        % plot
        plot(t, s1,'r','linewidth',1);
        hold on
        plot(t,s2,'b','linewidth',1);
        plot(t,s1+s2,'k','linewidth',1);
        % finalise axis
        axis([0 1 -2.2 2.2]);
        set(gca,'xtick',[],'ytick',[]);
    end
end
% label axis group
p(1,1).xlabel('time(unitless)');
p(1,1).ylabel('example data serier');
%% data set 2
source = 'XYZXYZ';
% an alternative way to access sub-panels is to first get a reference to
% the parent..
q = p(2);
% loop
for m = 1:6
    for n = 1:2
        % select axis -these two lines do the same thing above
        q(m,n).select();
        % prepare sample data
        data = randn(100,1)*0.4;
        % do stats
        stats = []; 
        stats.source = source(m);
        stats.binrange = [-1 1];
        stats.xtick = [-0.8:0.4:0.8];
        stats.ytick = [0 20];
        stats.bincens = -0.9:0.2:0.9;
        stats.values = data;
        stats.freq = hist(data,stats.bincens);
        stats.percfreq = stats.freq/length(data)*100;
        stats.percpeak = 30;
        % plot
        demopanel_minihist(stats,m==6,n==1);
    end
end
% label axis group
q.xlabel('data value (furlongs per fortnight');
q.ylabel('normalized frequency (%)');
% data set 3
p(1,2).select();
% prepare sample data
r1 = rand(100,1);
r2 = randn(100,1);
% plot 
plot(r1,r1+0.2*r2,'k.');
hold on
plot([0 1],[0 1],'r-');
legend('actual measurements','our predictions');
% finalise axis
xlabel('r_1');
ylabel('r_2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanel6 Packing is very flexible - it doesn't just do grids
clf 
% a) pack a pair of columns
% b) pack a bit into one of them and then pack some more
% c) pack into the other using absolute packing mode
% d) call select('all'), to show off the result
% a) create the root panel, and pack two columns. to pack columns instead
% of rows, we just pass "h" (horizontal) to pack() .
p = panel();
p.pack('h',2);
%p.select('all');
% b) pack some stuff into the left column
p(1).pack({1/6 1/6 1/6});
%p.select('all');

p(1).pack();
p(1).pack(); % p(1).pack({1/6 1/6 1/6 [] []});
% c) in other column, we will show how to do absolute mode packing, you can
% even place the child panel outside of its parent's area. just pass a
% 4-element row vector of [left bottom width height] to do absolute mode
% packing
p(2).pack({[-0.3 -0.01 1 0.4]})
%p(2).pack();
%p(2).pack({[0.2 0.61 0.6 0.4]});

% d) 
% use selectAll to quickly show the layout you've achieved.
% this commits all uncommitted panels as axis panels, so
% they can't be parents anymore (i.e. they can't have more
% children pack()ed into them).
p.select('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanel7 Panel gives figure-wide control over text properties
clf 
% a) Create a grid of panels
% b) Change some text properties
%% a)
% create a grid
p = panel();
p.pack(2,2);

% select all
p.select('all');

%% b) 
% if we set the properties on the root panel, they affect all its children
% and grandchildren.
p.fontname = 'Courier New';
p.fontsize = 10;
p.fontweight = 'normal'; % this is the default, anyway
% however, any child can override them, and the changes affect just that
% child(and its descendants).
p(2,2).fontsize = 14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanel8 Panels can be repacked from the command line
clf 
% a) Create a grid of panels, and show something in them
% b) Repack some of them, as if at the command line.
% a) 
% Create a 2x2 grid in gcf.
p = panel();
p.pack(2,2);
% have a look at p- all the child panels are currently uncommitted
p
% commit all the uncommitted panels as axis panels
p.select('all');
% b) During developement of a layout, you might find repack() useful
% repack one of the rows in the root panel
p(1).repack(0.3);
p.identify();
% repack one of the columns in one of the rows
p(1,1).repack(0.3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanel9 Panels can build complex layouts rapidly (HINTS on Margins!).
clf
% a) Build the layout from demopanel1, with annotation
% b) Add the content, so we can see what we are aiming for
% c) Show labelling of axis groups
% d) Add appropriate margins for this layout
%% a) 
% create panel
p = panel();
% let's start with two columns, one third and two thirds
p.pack('h',{1/3 2/3});
% then let's pack two rows into the first column, with the top row pretty
% big so we've room for some sub-panels
p(1).pack({2/3 []});
% now let's pack in those sub-panels
p(1,1).pack(3,2);
% finally, let's pack a grid of sub-panels into the right hand side too
p(2).pack(6,2);
%% b) % now let's populate those panels with axes full of data
% data set 1
for m = 1:3
    for n = 1:2
        % prepare sample data
        t = (0:99)/100;
        s1 = sin(t*2*pi*m);
        s2 = sin(t*2*pi*n*2);
        % select axis
        p(1,1,m,n).select();
        % q = p(1,1);
        %q(m,n).select();
        % plot
        plot(t,s1,'r','linewidth',1);
        hold on
        plot(t,s2,'b','linewidth',1);
        plot(t,s1+s2,'k','linewidth',1);
        % finalise axis
        axis([0 1 -2.2 2.2]);
        set(gca, 'xtick',[],'ytick',[]);
    end
end
% data set 2
source = 'XYZXYZ';
for m = 1:6
    for n = 1:2
        % select axis
        p(2,m,n).select();
        % prepare sample data
        data = randn(100,1)*0.4;
        % do stats
        stats = [];
        stats.source = source(m);
        stats.binrange = [-1 1];
        stats.xtick = [-0.8:0.4:0.8];
        stats.ytick = [0 20 40];
        stats.bincens = -0.9:0.2:0.9;
        stats.values = data;
        stats.freq = hist(data, stats.bincens);
        stats.percfreq = stats.freq/length(data)*100;
        stats.percpeak = 30;
        % plot
        demopanel_minihist(stats,m==6,n==1);
    end
end
% data set 3
p(1,2).select();
% prepare sample data
r1 = rand(100,1);
r2 = randn(100,1);
% plot
plot(r1,r1+0.2*r2,'k.')
hold on
plot([0 1],[0 1],'r-')
% finalise axis
xlabel('r_1')
ylabel('r_2')
legend('actual measurements','our predictions');
%% c) we can label parent panels (or "axis groups") just like labelling axis panels,
% except we have to use the method from panel, rather than the matlab call
% xlabel()
% label axis group
p(1,1).xlabel('time(unitless)');
p(1,1).ylabel('example data serier');
% we can also get a handle back to the label object, so that we can access
% its properties
% label axis group
h = p(2).xlabel('data value (furlongs per fortnight)');
p(2).ylabel('normalised frequency (%)')
% access propeties
get(h)

%% d) Let's do better with the margins
disp('These are the default margins-press any key to continue ...');
pause
%%% STEP 1: Tight internal margins
% tighten up all internal margins to the smallest margin
% we'll use anywhere (between the un-labelled sub-grids). this is usually a
% good starting point for any layout.
p.de.margin = 2;
% notice that we set the margin of all descendants of p, but the margin of
% p is not changed (p.de does not include p itself=, se there is still a
% margin from teh root panel p, to the figure edge. we can display this
% value:
disp(sprintf('p.margin is [%i %i %i %i ]',p.margin));
% the set p.fa(family) does include p, so p.fa is equal to {p.de and p}. if
% you could use p.fa.margin = 2 
% help panel/family panel/descendants
% pause
disp('We''ve tightened internal margins -press any key to continue ...');
pause

%%%% Step 2: Increase internal margins as required
% now let's space out the places we wnat to spaced out.
% remember that you can use p.identify() to get a nice indication of how to
% reference individual panels.
p(1,1).marginbottom = 12;
p(2).marginleft = 20;
% pause 
disp('We''ve increased two internal margins - press any key to continue...')
pause
%%%% STEP3 : Finalise margins with figure edges
% finally, let's sail as close to the wind as we dare for the final
% product, by trimming the root margin to the bone. eliminating any wasted
% whitespace like this is particularly helpful in exported image files
p.margin = [13 10 2 2];
% and let's set the global font properties, also,we can do this at any
% point, it doesn't have to be here
p.fontsize = 8;
% report 
% report
disp('We''ve now adjusted the figure edge margins (and reduced the fontsize), so we''re done.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanelA Panel builds immage files, not just on-screen figures 
clf
% a) Use demopanel1 to create a layout
% b) Export the result to file
% c) Export to different physical sizes
% d) Export at high quality
% e) Adjust margins
% f) Export using smoothing
% g) Export to EPS, rather than PNG
%% a)
% delegate
demopanel1
% see "help panel/export" for the full range of options
%% b) 
% the default sizing model for export targets a piece of paper. the default
% paper is A4, single column, with 20mm margins. the default aspect ratio
% is the golden ratio (landscape), therefore, if you provide only a
% filename you get this...
% do a default export
p.export('export_b'); % default resolution is 150DPI, so teh resulting file will look a bit scraggy,
% but it's a nice small file that is probably fine for laying out your
% document. note that we did not supply a file extension, so PNG format is
% assumed
%% c)
% one thing you might want to vary from figure to figure is the aspect
% ratio, the default (golden ratio) is a little short, here, so let's make
% it a touch taller.
p.export('export_c','-a1.4');
% the other thing is the column layout. we've exported to a
% single column, above - let's target a single column of a
% two-column layout.
p.export('export_c_c2', '-a1.4', '-c2');
% ach... that's never going to work, it doesn't fit in one
% column does it. this figure will have to span two columns,
% so let's leave it how it was.

% NB: here, we have used the "paper sizing model". if you
% prefer, you can use the "direct sizing model" and just
% specify width and height directly. see "help
% panel/export".
%% d) 
% when you're done drafting your document, you can bring up the export
% resolution to get a nice looking figure "-rp" means "publication
% resolution" (600DPI)
p.export('export_d','-a1.4','-rp');
%% e) Once exported at final resolution, let's put them as tight as we dare
% to reduce the whitespace
p.de.margin = 1;
p(1,1).marginbottom = 9;
p(2).marginleft = 12;
p.margin = [10 8 0.5 0.5];
p.export('export_e','-a1.4','-rp');
% NB: when the margins are this tight and the output
% resolution this high, you may notice small differences in
% layout between the on-screen renderer, the PNG renderer,
% and the EPS renderer.
%% f) that's now exported at 600DPI, which is fine for most purposes. however, the matlab renderer you
% use may not do nice anti-aliasing like some renderers. one way to
% mitigate this is to export at a higher DPI, btu that makes for a very
% large figure file. an alternative is to aks panel to render at a higher
% DPI but then to smotth back down to the specified resolution. you will
% have to wait a few seconds for the result, since rendering at these sizes
% takes time, here we'll smooth by factor 2. factor 4 takes even longer.
disp('rendering with smoothing, this may take longer ...');
p.export('export_f','-a1.4','-rp/2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanelB Panel can incorporate an existing axis
clf
% a) Create the root panel
% b) Create an axis yourself
% c) Pack an automatically created axis, and your own axis into the root
% panel
% a) create a column-pair layout, with 95% of the space given to the left
% hand panel
p = panel();
p.pack('h',{95 []});
% and put an axis in the left panel
h_axis = p(1).select();
% and  an image too
[X,Y,Z] = peaks(50);
surfc(X,Y,Z);
% b) sometimes you'll want to use some other function than panels to create one or more axes
% e.g. colorbar
h_colorbar_axis = colorbar('peer',h_axis);
% c) panel can manage the layout of these too
p(2).select(h_colorbar_axis);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demopanelC Recovering a Panel from a figure
clf
% a) Create a grid of panels and show something in them
% b) Recover the root panel from the figuer
% a) create a 2x2 grid in gcf
clf
p = panel();
p.pack(2,2);
% show dummy content
p.select('data');