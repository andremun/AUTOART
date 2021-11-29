%% Setting up the environment
global IMGTC IMGBIN IMGIND PRIM PATT Pr_PRIM I_PRIM

artworkfcn;
datadir = './data/';
imagedir = [datadir 'images/'];
nimgrow = 18; % Number of images per row
nimgcol = 17; % Number of images per col
nrowfig = 520; % Number of cols per image
ncolfig = 590; % Number of rows per image
nfigs = 306;

IMGTC = getfromfile([datadir 'raw_image_data.mat'],'IMGTC');
IMGIND = getfromfile([datadir 'raw_image_data.mat'],'IMGIND');
IMGBIN = getfromfile([datadir 'raw_image_data.mat'],'IMGBIN');
PRIM = getfromfile([datadir 'raw_image_data.mat'],'PRIM');
PATT = getfromfile([datadir 'raw_image_data.mat'],'PATT');
Pr_PRIM = getfromfile([datadir 'raw_image_data.mat'],'Pr_PRIM');
I_PRIM = getfromfile([datadir 'raw_image_data.mat'],'I_PRIM');

%% Colect the data from the Artwork. Produce the necessary images
Jnegtri = getfromfile([datadir 'result_triptych.mat'], 'Jnegtri');
Inegtri =  getfromfile([datadir 'result_triptych.mat'], 'Inegtri');

X = reshape(Inegtri(:,2),nimgrow,nimgcol);

clf;
imshow(rendercolor(X)); axis off;
print(gcf,'-dpng','-r300',[imagedir 'centre.png']);

clf;
subplot(2,2,1); imshow(rendercolor(X)); axis([0 3.*ncolfig 0 3.*nrowfig]); axis off;
subplot(2,2,2); imshow(renderindexed(X)); axis([0 3.*ncolfig 0 3.*nrowfig]); axis off;
subplot(2,2,3); imshow(renderbinary(X,true)); axis([0 3.*ncolfig 0 3.*nrowfig]); axis off;
subplot(2,2,4); imshow(renderbinary(I_PRIM(X),false)); axis([0 3.*ncolfig 0 3.*nrowfig]); axis off;
print(gcf,'-dpng','-r300',[imagedir 'representation.png']);

close all;
tiledlayout(1,3,'TileSpacing','compact');
for ii=1:3
    nexttile;
    imshow(rendercolor(reshape(Inegtri(:,ii),nimgrow,nimgcol)));
    axis off;
end
set(gcf,'WindowState','fullscreen');
print(gcf,'-dpng','-r300',[imagedir 'triptych.png']);

%% Collect some auxiliary images, such as the primitives and some patterns
close all;
for ii=1:size(PRIM,3)-2
    subplot(4,6,ii); imshow(PRIM(:,:,ii+2),'Border','tight'); axis off;
    title(['(' num2str(ii+2) ')']);
end
print(gcf,'-dpng','-r600 ',[imagedir 'primitives.png']);

clf;
idx = [2 3 12 21];
lbl = {'(A)','(B)','(C)','(D)'};
for ii=1:4
    subplot(2,2,ii);
    imshow([PRIM(:,:,PATT(1,1,idx(ii))) PRIM(:,:,PATT(1,2,idx(ii)));
            PRIM(:,:,PATT(2,1,idx(ii))) PRIM(:,:,PATT(2,2,idx(ii)))],'Border','tight');
    axis off;
    title(lbl(ii));
end
print(gcf,'-dpng','-r600',[imagedir 'patterns.png']);

%% Colect the data from the random generations. Produce the necessary images
if exist([datadir 'result_randart.mat'],'file')
    Jrand = getfromfile([datadir 'result_randart.mat'],'J');
    Irand = getfromfile([datadir 'result_randart.mat'],'I');
else
    Jrand = NaN.*ones(6,1e6);
    rng('default');
    Irand = zeros(nfigs,1e6);
    for ii=1:10 %
        Irand(:,ii) = randperm(nfigs);
    end
    a = 1:1000:1e6;
    b = a(2)-1:1000:1e6;
    for ii=1:1000
        range = a(ii):b(ii);
        Jrand(:,range) = randart([datadir 'rndresults/'],Irand(:,range),Jrand(:,range));
    end
    save([datadir 'result_randart.mat'],'J','I');
end

% [~,Imin] = min(Jrand,[],2);
% [~,Imax] = max(Jrand,[],2);
% 
% for ii=[2 5 6]
%     clf;
%     X = reshape(Irand(:,Imin(ii)),nimgrow,nimgcol);
%     subplot(1,2,1); imshow(rendercolor(X)); axis off;
%     
%     X = reshape(Irand(:,Imax(ii)),nimgrow,nimgcol);
%     subplot(1,2,2); imshow(rendercolor(X)); axis off;
%     
%     print(gcf,'-dpng','-r300',[imagedir 'extremes_randart_J' num2str(ii) '.png']);
% end

%% Collect the data from the optimized generations. Produce the necessary images
if exist([datadir 'result_autoart.mat'],'file')
    Jopt = getfromfile([datadir 'result_autoart.mat'],'Jopt');
    Iopt = getfromfile([datadir 'result_autoart.mat'],'Iopt');
    Mopt = getfromfile([datadir 'result_autoart.mat'],'Mopt');
else
    Jopt = NaN.*ones(1,60);
    Iopt = zeros(nfigs,60);
    Mopt = zeros(3,60);
    Ftype = [2 5 6];
    acc = 0;
    for ii=1:2
        for jj=1:length(Ftype)
            for kk=1:10
                acc = acc+1;
                Mopt(:,acc) = [kk Ftype(jj) ii==2];
                filename = [datadir 'autoresults/result_S' num2str(kk) ...
                                                      '_E' num2str(Ftype(jj)) ...
                                                      '_M' num2str(ii==2) '.mat'];
                % disp(filename)
                if ~exist(filename,'file')
                    continue;
                end
                X = getfromfile(filename,'X');
                Iopt(:,acc) = X(:);
                J = getfromfile(filename,'J');
                Jopt(acc) = J(end);
            end
        end
    end
    
end

%%
Jrand = Jrand([2 5 6],:);
Jnegtri = Jnegtri([2 5 6],:);
Jopt = [Jopt([1:10 31:40]);
        Jopt([11:20 41:50]);
        Jopt([21:30 51:60])];
Jmax = max([Jrand Jopt Jnegtri],[],2);
Jmin = min([Jrand Jopt Jnegtri],[],2);
Jrandn = (Jrand-Jmin)./(Jmax-Jmin);
Jnegtrin = (Jnegtri-Jmin)./(Jmax-Jmin);
Joptn = (Jopt-Jmin)./(Jmax-Jmin);

%%
[~,idxmin] = min(Jopt,[],2);
idxmin = idxmin'+[0 10 20];
[~,idxmax] = max(Jopt,[],2);
idxmax = idxmax'+[20 30 40];
IMGTC_BACKUP = IMGTC;

for ii=1:3
    clf;
    filename = [datadir 'autoresults/result_S' num2str(Mopt(1,idxmin(ii))) ...
                                          '_E' num2str(Mopt(2,idxmin(ii))) ...
                                          '_M' num2str(Mopt(3,idxmin(ii))) '.mat'];
    IMGTC = getfromfile(filename,'IMGTC'); %#ok<NASGU>
    X = reshape(Iopt(:,idxmin(ii)),nimgrow,nimgcol);
    tiledlayout(1,3,'TileSpacing','compact');
    
    nexttile; imshow(rendercolor(X)); axis off;
    
    IMGTC = IMGTC_BACKUP; %#ok<NASGU>
    nexttile; imshow(rendercolor(reshape(Inegtri(:,2),nimgrow,nimgcol))); axis off;
    
    
    filename = [datadir 'autoresults/result_S' num2str(Mopt(1,idxmax(ii))) ...
                                          '_E' num2str(Mopt(2,idxmax(ii))) ...
                                          '_M' num2str(Mopt(3,idxmax(ii))) '.mat'];
    IMGTC = getfromfile(filename,'IMGTC');
    X = reshape(Iopt(:,idxmax(ii)),nimgrow,nimgcol);
    nexttile; imshow(rendercolor(X)); axis off;
    
    set(gcf,'WindowState','fullscreen');
    print(gcf,'-dpng','-r300',[imagedir 'extremes_autoart_J' num2str(ii) '.png']);
end

%%
close all;
h = zeros(4,1);
labels = {'MI','CP','PI'};
for ii=1:3
    subplot(3,1,ii);
    h(1) = histogram(Jrandn(ii,:),50,'Normalization','probability');
    h(2) = line([Jnegtrin(ii,1) Jnegtrin(ii,1)],[0 2e5], ...
                                          'LineStyle', '--', ...
                                          'Color', 'r');
    h(3) = line([Jnegtrin(ii,2) Jnegtrin(ii,2)],[0 2e5], ...
                                        'LineStyle', '--', ...
                                        'Color', 'k');
    h(4) = line([Jnegtrin(ii,3) Jnegtrin(ii,3)],[0 2e5], ...
                                        'LineStyle', '--', ...
                                        'Color', 'g');
    h(6) = line(Joptn(ii,1:10),0.15.*ones(10,1), 'LineStyle', 'none', ...
                                                 'Marker', 'x' ,...
                                                 'Color', 'm');
    h(5) = line(Joptn(ii,11:20),0.15.*ones(10,1), 'LineStyle', 'none', ...
                                                  'Marker', 'x' ,...
                                                  'Color', 'b');
    ylabel(['Pr_{' labels{ii} '}']);
    legend(h,{'RAND','LEFT','CENTRE','RIGHT','OPT_{max}','OPT_{min}'},'location','NorthEastOutside');
    axis([-0.05 1.05 0 0.20]);
end
print(gcf,'-dpng',[imagedir 'histograms.png']);


