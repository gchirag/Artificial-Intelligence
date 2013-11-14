function [Y, R, E, M, N] = Isomap(images, ks, labels, D, n_fcn, n_size, options); 

% ISOMAP   Computes Isomap embedding using the algorithm of 
%             Tenenbaum, de Silva, and Langford (2000). 
%
% [Y, R, E] = isomap(D, n_fcn, n_size, options); 
%
% Input:
%    D = N x N matrix of distances (where N is the number of data points)
%    n_fcn = neighborhood function ('epsilon' or 'k') 
%    n_size = neighborhood size (value for epsilon or k) 
%
%    options.dims = (row) vector of embedding dimensionalities to use
%                        (1:10 = default)
%    options.comp = which connected component to embed, if more than one. 
%                        (1 = largest (default), 2 = second largest, ...)
%    options.display = plot residual variance and 2-D embedding?
%                        (1 = yes (default), 0 = no)
%    options.overlay = overlay graph on 2-D embedding?  
%                        (1 = yes (default), 0 = no)
%    options.verbose = display progress reports? 
%                        (1 = yes (default), 0 = no)
%
% Output: 
%    Y = Y.coords is a cell array, with coordinates for d-dimensional embeddings
%         in Y.coords{d}.  Y.index contains the indices of the points embedded.
%    R = residual variances for embeddings in Y
%    E = edge matrix for neighborhood graph
%

%    BEGIN COPYRIGHT NOTICE
%
%    Isomap code -- (c) 1998-2000 Josh Tenenbaum
%
%    This code is provided as is, with no guarantees except that 
%    bugs are almost surely present.  Published reports of research 
%    using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%      J. B. Tenenbaum, V. de Silva, J. C. Langford (2000).  A global
%      geometric framework for nonlinear dimensionality reduction.  
%      Science 290 (5500): 2319-2323, 22 December 2000.  
%
%    Comments and bug reports are welcome.  Email to jbt@psych.stanford.edu. 
%    I would also appreciate hearing about how you used this code, 
%    improvements that you have made to it, or translations into other
%    languages.    
%
%    You are free to modify, extend or distribute this code, as long 
%    as this copyright notice is included whole and unchanged.  
%
%    END COPYRIGHT NOTICE


%%%%% Step 0: Initialization and Parameters %%%%%

N = size(D,1); 
if ~(N==size(D,2))
     error('D must be a square matrix'); 
end; 
if n_fcn=='k'
     K = n_size; 
     if ~(K==round(K))
         error('Number of neighbors for k method must be an integer');
     end
elseif n_fcn=='epsilon'
     epsilon = n_size; 
else 
     error('Neighborhood function must be either epsilon or k'); 
end
if nargin < 6
     error('Too few input arguments'); 
elseif nargin < 7
     options = struct('dims',1:10,'overlay',1,'comp',1,'display',1,'verbose',1); 
end
INF =  1000*max(max(D))*N;  %% effectively infinite distance

if ~isfield(options,'dims')
     options.dims = 1:10; 
end
if ~isfield(options,'overlay')
     options.overlay = 1; 
end
if ~isfield(options,'comp')
     options.comp = 1; 
end
if ~isfield(options,'display')
     options.display = 1; 
end
if ~isfield(options,'verbose')
     options.verbose = 1; 
end
dims = options.dims; 
comp = options.comp; 
overlay = options.overlay; 
displ = options.display; 
verbose = options.verbose; 

Y.coords = cell(length(dims),1); 
R = zeros(1,length(dims)); 

%%%%% Step 1: Construct neighborhood graph %%%%%
disp('Constructing neighborhood graph...'); 

if n_fcn == 'k'
     [tmp, ind] = sort(D); 
     for i=1:N
          D(i,ind((2+K):end,i)) = INF; 
     end
elseif n_fcn == 'epsilon'
     warning off    %% Next line causes an unnecessary warning, so turn it off
     D =  D./(D<=epsilon); 
     D = min(D,INF); 
     warning on
end

D = min(D,D');    %% Make sure distance matrix is symmetric

if (overlay == 1)
     E = int8(1-(D==INF));  %%  Edge information for subsequent graph overlay
end

% Finite entries in D now correspond to distances between neighboring points. 
% Infinite entries (really, equal to INF) in D now correspond to 
%   non-neighoring points. 

%%%%% Step 2: Compute shortest paths %%%%%
disp('Computing shortest paths...'); 

% We use Floyd's algorithm, which produces the best performance in Matlab. 
% Dijkstra's algorithm is significantly more efficient for sparse graphs, 
% but requires for-loops that are very slow to run in Matlab.  A significantly 
% faster implementation of Isomap that calls a MEX file for Dijkstra's 
% algorithm can be found in isomap2.m (and the accompanying files
% dijkstra.c and dijkstra.dll). 

tic; 
for k=1:N
     D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
     if ((verbose == 1) & (rem(k,20) == 0)) 
          disp([' Iteration: ', num2str(k), '     Estimated time to completion: ', num2str((N-k)*toc/k/60), ' minutes']); 
     end
end

%%%%% Remove outliers from graph %%%%%
disp('Checking for outliers...'); 
n_connect = sum(~(D==INF));        %% number of points each point connects to
[tmp, firsts] = min(D==INF);       %% first point each point connects to
[comps, I, J] = unique(firsts);    %% represent each connected component once
size_comps = n_connect(comps);     %% size of each connected component
[tmp, comp_order] = sort(size_comps);  %% sort connected components by size
comps = comps(comp_order(end:-1:1));    
size_comps = size_comps(comp_order(end:-1:1)); 
n_comps = length(comps);               %% number of connected components
if (comp>n_comps)                
     comp=1;                              %% default: use largest component
end
disp(['  Number of connected components in graph: ', num2str(n_comps)]); 
disp(['  Embedding component ', num2str(comp), ' with ', num2str(size_comps(comp)), ' points.']); 
Y.index = find(firsts==comps(comp)); 

D = D(Y.index, Y.index); 
N = length(Y.index); 

%%%%% Step 3: Construct low-dimensional embeddings (Classical MDS) %%%%%
disp('Constructing low-dimensional embeddings (Classical MDS)...'); 

opt.disp = 0; 
[vec, val] = eigs(-.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2)), max(dims), 'LR', opt); 

h = real(diag(val)); 
[foo,sorth] = sort(h);  sorth = sorth(end:-1:1); 
val = real(diag(val(sorth,sorth))); 
vec = vec(:,sorth); 

D = reshape(D,N^2,1); 
for di = 1:length(dims)
     if (dims(di)<=N)
         Y.coords{di} = real(vec(:,1:dims(di)).*(ones(N,1)*sqrt(val(1:dims(di)))'))'; 
         r2 = 1-corrcoef(reshape(real(L2_distance(Y.coords{di}, Y.coords{di})),N^2,1),D).^2; 
         R(di) = r2(2,1); 
         if (verbose == 1)
             disp(['  Isomap on ', num2str(N), ' points with dimensionality ', num2str(dims(di)), '  --> residual variance = ', num2str(R(di))]); 
         end
     end
end

clear D; 

%%%%%%%%%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%

if (displ==1)

     %%%%% Plot two-dimensional configuration %%%%%
     twod = find(dims==2); 
     if ~isempty(twod)
         fig=figure;
         hold on;
         M=Y.coords{twod}(1,:);
         N=Y.coords{twod}(2,:);
         poses(1,:)=M;
		 poses(2,:)=N;
       %  plot(Y.coords{twod}(1,:), Y.coords{twod}(2,:), 'ro'); 
			plot(0,0,'yo');
			plot(0,0,'m^');
			plot(0,0,'c^');
			plot(0,0,'ro');
			plot(0,0,'go');
			plot(0,0,'ko');
			plot(0,0,'r*');
			plot(0,0,'g*');
			plot(0,0,'c+');
			plot(0,0,'b*');         for j = 1: size(Y.coords{twod} , 2);
            if ( labels(j) == 0)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'yo');
            end

            if ( labels(j) == 1)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'm^');
            end
            if ( labels(j) == 2)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'c^');
            end

            if ( labels(j) == 3)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'ro');
            end
            if ( labels(j) == 4)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'go');
            end

            if ( labels(j) == 5)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'ko');
            end
            if ( labels(j) == 6)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'r*');
            end

            if ( labels(j) == 7)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'g*');
            end
            if ( labels(j) == 8)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'c+');
            end

            if ( labels(j) == 9)
                plot(Y.coords{twod}(1,j) , Y.coords{twod}(2,j) , 'b*');
            end
         end
         legend('0','1','2','3','4','5','6','7','8','9');
         hold on
         plot(poses(1,ks),poses(2,ks),'bo');
         hold on
         scale1 = range(poses(1,:))/800;
         scale2 = range(poses(2,:))/1600;

		x=zeros(28,28);

	for p=1:size(ks,2)

    k=ks(p);

    for i=1:28

        x(i,:)=images((i-1)*28+1:i*28,k);

    end
    colormap(gray);

    xc=poses(1,k);

    yc=poses(2,k);

  rotate(image([xc xc-28*2*scale1],[yc yc-28*2*scale2], x),[0,0,1],90,[xc,yc,0]);
  hold on
end
         if (overlay == 1)
             %gplot(E(Y.index, Y.index), [Y.coords{twod}(1,:); Y.coords{twod}(2,:)]'); 
             title('Two-dimensional Isomap embedding (with neighborhood graph).'); 
         else
             title('Two-dimensional Isomap.'); 
         end
         hold off;
 		 saveas(fig,'Euclidean_all_digits.jpg');

     end
end

