function [V,C,XY]=VoronoiLimit(varargin)
% --------------------------------------------------------------
% [V,C,XY]=VoronoiLimit(x,y,additional_variables)
% Provides the Voronoi decomposition of a set of (x,y) data, but with all
% vertices limited to the boundary created by the data itself.
% V contains all vertices and C contains all vertices for each individual
% point. That is: V(C{ij},:) will give you the vertices of the ij'th data
% point. The order of polygon vertices are given in a counter-clockwise
% manner. XY contains updated xy coordinates as limited by any input boundaries.
%
% Addition variables:
% 'bs_ext':  Describe an arbitrary external boundary by giving an xy matrix of size (n,1) where n are number of vertices.
% 'bs_int':  Describe any number of arbitrary internal boundaries by giving a cell structure of M xy matrices of size
%            (Ni,1) where M are number of internal boundaries and Ni are number of vertices in the respective boundaries.
% 'figure':  output figure ('on'/'off'. Default='on').
%
% Run with no input to see example.
% Check the code for this case to see examples of the use of external and internal boundaries.
%
% Requires the Polybool function of the mapping toolbox to run!.
% I recommend the tool 'export_fig' for exporting figures. It can be found here:
% http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
%
% Made by: Jakob Sievers, Jakob.Sievers@gmail.com
% Edits:
% 25 March 2015: Revised to support internal and external boundaries
% 29 March 2015: Bugfixes
% 31 March 2015: Bugfixes
% 26  June 2015: Bugfixes + addition of a triangle external boundary example 
% 29  June 2015: Output vertex order set to counter-clockwise
% 7   July 2015: Treatment of xy input points around boundaries improved (See example figure)
% 21  July 2015: Routine improved to allow for proper handling of cells which are split into two, or more, closed polygons. In these cases the original cell is split into an appropriate number of new cells. 
% --------------------------------------------------------------
warning('off','map:polygon:noExternalContours');

   
if nargin==0
    val=600;
    x=rand(val,1);
    y=rand(val,1);
    XY=unique([x,y],'rows');
    x=XY(:,1);
    y=XY(:,2);
    
    %EXTERNAL BOUNDARIES
    ButtonName = questdlg('Choose external boundary example:','','Irregular pentagon', 'Triangle', 'Irregular pentagon');
    switch ButtonName,
        case 'Irregular pentagon',
            bs_ext=[min(x)-std(x)/2 min(x)-std(x)/2 0.65 max(x)+std(x)/2 max(x)+std(x)/2 min(x)-std(x)/2;min(y)-std(y)/2 max(y)+std(y)/2 max(y)+std(y)/2 .65 min(y)-std(y)/2 min(y)-std(y)/2]';
        case 'Triangle',
            bs_ext=[-.8 .5 1.80 -.8;-.05 1.7 -.05 -.05]';
            save bs_ext
    end
    save bs_ext
    %INTERNAL OBJECTS
    bs_int=cell(3,1);
    rat=1.5;
    % rectangle
    bs_int{1}=[min(x)+(std(x)*rat) min(x)+(std(x)*rat) max(x)-std(x) max(x)-std(x) min(x)+(std(x)*rat);min(y)+std(y) max(y)-std(y) max(y)-std(y) min(y)+std(y) min(y)+std(y)]';
    t = linspace(0,2*pi)';
    % circle 1
    xc=.25;
    yc=.7;
    rad=.10;
    bs_int{2}=[(cos(t)*rad)+xc (sin(t)*rad)+yc];
    % circle 2
    xc=.4;
    yc=.3;
    rad=.16;
    bs_int{3}=[(cos(t)*rad)+xc (sin(t)*rad)+yc];
    fig='on';
else
    x=varargin{1}(:);
    y=varargin{2}(:);
    XY=unique([x,y],'rows');
    x=XY(:,1);
    y=XY(:,2);
    for ii=3:2:nargin
        if strcmp(varargin{ii},'bs_ext')
            bs_ext=varargin{ii+1};
        elseif strcmp(varargin{ii},'bs_int')
            bs_int=varargin{ii+1};
        elseif strcmp(varargin{ii},'figure')
            fig=varargin{ii+1};
        end
    end
    if exist('fig','var')==0
        fig='on';
    end
end


x=x(:);
y=y(:);
rx=[min(x) max(x)];
ry=[min(y) max(y)];

bnd=[rx ry]; %data bounds

crs=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]); %data boundary corners

if exist('bs_ext','var')
    crs=bs_ext;
end
% if exist('bs_int','var')
%     for ii=1:length(bs_int)
%         inpol=inpolygon(x,y,bs_int{ii}(:,1),bs_int{ii}(:,2));
%         x(inpol)=[];
%         y(inpol)=[];
%     end
% end

if ~any(size(x)==1) || ~any(size(y)==1) || numel(x)==1 || numel(y)==1
    disp('Input vectors should be single rows or columns')
    return
end


dt=delaunayTriangulation(x(:),y(:));
[V,C]=voronoiDiagram(dt); %This structure gives vertices for each individual point but is missing all "infinite" vertices
[vx,vy]=voronoi(x,y); %This structure includes the "infinite" vertices but provides everything as a completele list of vertices rather than individually for each point.
%Hence we need to add the missing vertices from vx and vy to the V and C structure.
vxyl=[vx(:) vy(:)];
xix=ones(size(vx));


%values provided by voronoiDiagram may be an infinitesimal fraction off
%relative to those provided by "voronoi". Hence we need to make sure all
%values in V are similar to those located in vxyl.
vals=unique(vxyl(:));
for ik=1:length(vals)
    ix=find(V(:)==vals(ik));
    if ~isempty(ix)
        V(ix)=vals(ik);
    end
end
lV0=length(V);

%Find missing points that should be added to existing V/C structure
for ii=1:length(vxyl)
    fix=find(V(:,1)==vxyl(ii,1));
    if ~isempty(fix)
        if any(V(fix,2)==vxyl(ii,2))
            xix(ii)=0;
        end
    end
end

mix=find(xix==1)./2; %index of missing values
lmix=length(mix);
mvx=vx(2,mix); %missing vx
mvy=vy(2,mix); %missing vy
mv=[mvx',mvy'];
cpx=vx(1,mix); %connector point x (connects between outer missing points and inner existing points in V/C)
cpy=vy(1,mix); %connector point y (connects between outer missing points and inner existing points in V/C)

ctr=0;
mv2=[];
cpVixt=cell(lmix,1); %connector points, index in V structure
for ii=1:lmix
    if any(V(:,1)==cpx(ii) & V(:,2)==cpy(ii))
        cpVixt{ii}=find(V(:,1)==cpx(ii) & V(:,2)==cpy(ii));
        lval=length(cpVixt{ii});
        if lval==1
            ctr=ctr+1;
            mv2(ctr,:)=mv(ii,:);
        elseif lval>1
            ctr=ctr+1;
            mv2(ctr:ctr+lval-1,:)=[ones(lval,1).*mv(ii,1) ones(lval,1).*mv(ii,2)];
            ctr=ctr+lval-1;
        end
    end
end
cpVixt=cell2mat(cpVixt);


V=[V;mv2]; %add points to V structure

%Addition-routine: addition of missing points (mvx,mvy) to individual vertice-polygons (C)
for ij=1:length(C)
    if any(C{ij}==1)
        ixa=find(cpVixt==C{ij}(2));
        ixb=find(cpVixt==C{ij}(end));
        if  length(C{ij})<3
            C{ij}(1)=lV0+ixa(1);
            C{ij}=[C{ij},lV0+ixa(2)];
        else
            if length(ixa)==1 && length(ixb)==1
                C{ij}(1)=lV0+ixa;
                C{ij}=[C{ij},lV0+ixb];
            elseif length(ixa)==2 && length(ixb)==1
                C{ij}=[C{ij},lV0+ixb];
                [~,minix]=min(sqrt((V(C{ij}(end),1)-V(lV0+ixa,1)).^2+(V(C{ij}(end),2)-V(lV0+ixa,2)).^2));
                C{ij}(1)=lV0+ixa(minix);
            elseif length(ixa)==1 && length(ixb)==2
                C{ij}(1)=lV0+ixa;
                [~,minix]=min(sqrt((V(C{ij}(1),1)-V(lV0+ixb,1)).^2+(V(C{ij}(1),2)-V(lV0+ixb,2)).^2));
                C{ij}=[C{ij},lV0+ixb(minix)];
            elseif length(ixa)==2 && length(ixb)==2
                dist1=sqrt((x(ij)-V(lV0+ixa,1)).^2+(y(ij)-V(lV0+ixa,2)).^2);
                dist2=sqrt((x(ij)-V(lV0+ixb,1)).^2+(y(ij)-V(lV0+ixb,2)).^2);
                if diff(dist1)==0 && diff(dist2)==0
                    minix1=1;
                    minix2=2;
                else
                    [~,minix1]=min(dist1);
                    [~,minix2]=min(dist2);
                end
                C{ij}(1)=lV0+ixa(minix1);
                C{ij}=[C{ij},lV0+ixb(minix2)];
            end
        end
    end
end

%Polybool for restriction of polygons to domain.
C1=C; %Do this analysis based on old vertice descriptions to avoid problems
allVixinp=inpolygon(V(:,1),V(:,2),crs(:,1),crs(:,2)); %determine which points in V that are within the data boundaries.
isemp=false(length(C),1);
for ij=1:length(C)
    if sum(allVixinp(C{ij}))~=length(C{ij})
        if license('test','MAP_Toolbox')
            [xb, yb] = polybool('intersection',crs(:,1),crs(:,2),V(C1{ij},1),V(C1{ij},2));
        else
            [xb, yb] = newPolybool(crs(:,1),crs(:,2),V(C1{ij},1),V(C1{ij},2),1);
        end
        ix=nan(1,length(xb));
        for il=1:length(xb)
            if any(V(:,1)==xb(il)) && any(V(:,2)==yb(il))
                ix1=find(V(:,1)==xb(il));
                ix2=find(V(:,2)==yb(il));
                for ib=1:length(ix1)
                    if any(ix1(ib)==ix2)
                        ix(il)=ix1(ib);
                    end
                end
                if isnan(ix(il))==1
                    lv=length(V);
                    V(lv+1,1)=xb(il);
                    V(lv+1,2)=yb(il);
                    allVixinp(lv+1)=1;
                    ix(il)=lv+1;
                end
            else
                lv=length(V);
                V(lv+1,1)=xb(il);
                V(lv+1,2)=yb(il);
                allVixinp(lv+1)=1;
                ix(il)=lv+1;
            end
        end
        C{ij}=ix;
    end
    if isempty(C{ij})
        isemp(ij)=true;
    end
end
if any(isemp)
    C(isemp)=[];
    XY(isemp,:)=[];
end


%adjust polygons to the presence of internal boundaries
if exist('bs_int','var')
    isemp=false(length(C),length(bs_int));
    for ii=1:length(bs_int)
        V2=nan(length(V)*10,2);
        C2=cell(length(C),1);
        ctr=1;
        for ij=1:length(C)
            if license('test','MAP_Toolbox')
                [pbx,pby]=polybool('subtraction',V(C{ij},1),V(C{ij},2),bs_int{ii}(:,1),bs_int{ii}(:,2));
            else
                [pbx,pby] = newPolybool(V(C{ij},1),V(C{ij},2),bs_int{ii}(:,1),bs_int{ii}(:,2),3);
            end
            if ~isempty(pbx)
                C2{ij}=(ctr:ctr+length(pbx)-1)';
                C2{ij}=[C2{ij} ones(size(C2{ij}))*ij];
                V2(ctr:ctr+length(pbx)-1,:)=[pbx pby];
                ctr=ctr+length(pbx);
            end
        end
        V=V2(1:ctr-1,:);
        C=C2;
        for ij=1:length(C)
            if isempty(C{ij})
                isemp(ij,ii)=true;
            else
                C{ij}=(C{ij}(:,1))';
            end
        end
    end
    if any(any(isemp'))
        C(any(isemp'))=[];
        XY(any(isemp'),:)=[];
    end
end


%remove spurious double-entires in C/V structure
epsx=eps(max(abs(V(isinf(V)==0))));
for ih=1:length(C)
    VC=V(C{ih},:);
    TMAT=true(size(VC,1));
    for ii=1:size(VC,1)
        for ij=1:size(VC,1)
            TMAT(ii,ij)=all(abs(VC(ii,:)-VC(ij,:))<=epsx);
        end
    end
    TMAT=TMAT-eye(size(TMAT));
    if any(TMAT(:)==1)
        if all(abs(V(C{ih}(1),:)-V(C{ih}(end),:))<=epsx)
            C{ih}(end)=[];
        end
        ctr=0;
        while ctr<length(C{ih})-1
            ctr=ctr+1;
            if all(abs(V(C{ih}(ctr),:)-V(C{ih}(ctr+1),:))<=epsx)
                C{ih}(ctr+1)=[];
            end
        end
    end
    C{ih}=C{ih}';
end


TMAT=cell(length(V)-1,1);
Vt=V;
idx1=(1:length(V))';
idx2=(1:length(V))';
for ii=1:length(V)-1
    Vt=[Vt(2:end,:);Vt(1,:)];
    idx2=[idx2(2:end);idx2(1)];
    TMATt=find(all(abs(V-Vt)<=epsx,2));
    TMAT{ii}=[idx1(TMATt) idx2(TMATt)];
end
TMATf=unique(sort(cell2mat(TMAT),2),'rows');
if ~isempty(TMATf)
    for ii=1:size(TMATf,1)
        for ij=1:length(C)
            C{ij}(C{ij}==TMATf(ii,2))=TMATf(ii,1);
        end
    end
end



%remove V-entries which are now unused by C
index_rem=true(size(V,1),1);
Ctot=unique(cell2mat(C));
index_rem(Ctot)=false;
index_rem=find(index_rem);
while ~isempty(index_rem)
    for ij=1:length(C)
        ixf=find(C{ij}>index_rem(1));
        if ~isempty(ixf)
            C{ij}(ixf)=C{ij}(ixf)-1;
        end
    end
    V(index_rem(1),:)=[];
    index_rem=true(size(V,1),1);
    Ctot=unique(cell2mat(C));
    index_rem(Ctot)=false;
    index_rem=find(index_rem);
end

%Check and repair cells that have been split into closed sub-cells by input boundaries
Csplit=cell(length(C),1);
XYsplit=cell(length(C),1);
splitlog=false(length(C),1);
for ij=1:length(C)
    [xClosed, yClosed] = closePolygonParts(V(C{ij},1),V(C{ij},2));
    if any(isnan(xClosed))
        splitlog(ij)=true;
        ix=find(~isnan(xClosed));
        diffix=diff(ix)>1;
        NUMcell=sum(isnan(xClosed))+1;
        Csplit{ij}=cell(NUMcell,1);
        XYsplit{ij}=nan(NUMcell,2);
        C_temp=C{ij};
        ix_begin=1;
        for ik=1:NUMcell 
            cs_diffix=cumsum(diffix);
            if ik>1
               ix_begin=2; 
            end
            ix_end=find(cs_diffix>0,1,'first');
            if isempty(ix_end)
               ix_end=length(xClosed); 
            end
            Csplit{ij}{ik}=C_temp(ix_begin:ix_end);
            inpol=inpolygon(XY(ij,1),XY(ij,2),xClosed(ix_begin:ix_end),yClosed(ix_begin:ix_end));
            if inpol==0
                XYsplit{ij}(ik,:)=[mean(xClosed(ix_begin:ix_end)) mean(yClosed(ix_begin:ix_end))];
            else
                XYsplit{ij}(ik,:)=XY(ij,:);
            end
            if ik<NUMcell
                C_temp(ix_begin:ix_end)=[];
                diffix(ix_begin:ix_end)=[];
                xClosed(ix_begin:ix_end)=[];
                yClosed(ix_begin:ix_end)=[];
            end
        end
    end
end
if any(splitlog)
    ix_splitlog=find(splitlog);
    ix_splitlog0=ix_splitlog;
    for ij=1:length(ix_splitlog)
        if ix_splitlog(ij)==1
            C=[Csplit{ix_splitlog(ij)};C(2:end)];
            XY=[XYsplit{ix_splitlog(ij)};XY(2:end,:)];
        elseif ix_splitlog(ij)==length(C)
            C=[C(1:end-1);Csplit{ix_splitlog(ij)}];
            XY=[XY(1:end-1,:);XYsplit{ix_splitlog(ij)}];
        else
            C=[C(1:ix_splitlog(ij)-1);Csplit{ix_splitlog0(ij)};C(ix_splitlog(ij)+1:end)];
            XY=[XY(1:ix_splitlog(ij)-1,:);XYsplit{ix_splitlog0(ij)};XY(ix_splitlog(ij)+1:end,:)];
            if ij<length(ix_splitlog)
                ix_splitlog(ij+1:end)=ix_splitlog(ij+1:end)+(length(Csplit{ix_splitlog0(ij)})-1);
            end
        end
    end
end

%ensure that all polygon vertex groups are given in counter-clockwise order
for ih=1:length(C)
   if ispolycw(V(C{ih},1),V(C{ih},2)) 
       C{ih}=flipud(C{ih});
   end
end

%close polygons for the purpose of plotting
C2=C;
for ih=1:length(C2)
    if C2{ih}(1)~=C2{ih}(end)
        C2{ih}=[C2{ih};C2{ih}(1)];
    end
end

% create and output figure
if exist('fig','var')
    if strcmp(fig,'on')
        figure
        set(gcf,'position',get(0,'screensize'),'color','w')
        set(gca,'box','on')
        hold on
        plot(x,y,'.k')
        if any(splitlog)
            for ij=1:length(ix_splitlog0)
                plot(XYsplit{ix_splitlog0(ij)}(:,1),XYsplit{ix_splitlog0(ij)}(:,2),'*r')
                plot(XYsplit{ix_splitlog0(ij)}(:,1),XYsplit{ix_splitlog0(ij)}(:,2),'or','markersize',8)
            end
        end
        voronoi(x,y)
        for id=1:length(C2)
            plot(V(C2{id},1),V(C2{id},2),'-r')
        end
        grid on
        axis tight
        axis square
        if nargin==0
            axis equal
        end
        ax=axis;
        dx=(ax(2)-ax(1))/10;
        dy=(ax(4)-ax(3))/10;
        axis([ax(1)-dx ax(2)+dx ax(3)-dy ax(4)+dy])
        title({'Original Voronoi Decomposition ({\color{blue}blue})';'New limited Voronoi Decomposition ({\color{red}red})'},'fontsize',16,'fontweight','bold')
        if exist('bs_int','var')
            for ii=1:length(bs_int)
                text(mean(unique(bs_int{ii}(:,1))),mean(unique(bs_int{ii}(:,2))),num2str(ii),'fontsize',30,'fontweight','bold','horizontalalignment','center')
            end
        end
    end
end
% % %     export_fig([pwd,'\VoronoiLimit_example.jpg'],'-r300');
% % 
% %     
% % end