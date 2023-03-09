function out = usaExtract(varargin)
% Check if the function has arguments

if(nargin==0)
    % Load set of coordinates
    load('../data/usa_coord.mat')
else
    usa = varargin{1};
end

% Create matrices of long and lat 
lat = repmat(flip(-89.5:89.5)',1,360);
long = repmat(-179.5:179.5,180,1);

% Loop over USA coordinates
N = size(usa,1);
X = [lat(:),long(:)];



% Extract vector
out = zeros(size(lat(:)));

% Loop getting 
for j= 1:N
    % Find indices that minimize
    [~,I] = min(sum((X-usa(j,:)).^2,2));
    % Replace on out vector as part of USA
    out(I) = true;
end

out = reshape(out,180,360);

end