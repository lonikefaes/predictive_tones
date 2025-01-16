function [gx, gy, gz] = gradfourier(potential),

[nx, ny, nz] = size(potential);
%% fourier
potential(isnan(potential)) = 0;
fourier = fftshift(fftn(potential)) .* tukeyWindow([nx ny nz]);
fourier = ifftshift(fourier);

hx = ceil(nx / 2) - 1;
ftdiff = (2i * pi / sqrt(nx)) * (0:hx);
ftdiff(nx:-1:nx - hx + 1) = -ftdiff(2 : hx + 1); % correct conjugate symmetry
fourierX = -fourier .* repmat(ftdiff', [1, ny, nz]);

hy = ceil(ny / 2) - 1;
ftdiff = (2i * pi / sqrt(ny)) * (0:hy);
ftdiff(ny:-1:ny - hy + 1) = -ftdiff(2 : hy + 1); % correct conjugate symmetry
fourierY = -fourier .* repmat(ftdiff,  [nx, 1, nz]);

hz = ceil(nz / 2) - 1;
ftdiff = (2i * pi / sqrt(nz)) * (0:hz);
ftdiff(nz:-1:nz - hz + 1) = -ftdiff(2 : hz + 1); % correct conjugate symmetry
fourierZ = -fourier .* permute(repmat(ftdiff', [1, nx, ny]), [2, 3, 1]);

gx = real(ifftn(fourierX));
gy = real(ifftn(fourierY));
gz = real(ifftn(fourierZ));
end

function window = tukeyWindow(windowSize)

r = 0.5;
numberOfDimensions = length(windowSize);
window1D = cell(numberOfDimensions, 1);
for i = 1:numberOfDimensions
    window1D{i} = tukeywin(windowSize(i) + mod(windowSize(i), 2), r);
    window1D{i} = window1D{i}(1:end - mod(windowSize(i), 2));
end
[window1D{:}] = meshgrid(window1D{:});
window = reshape([window1D{:}], [windowSize(2), windowSize(1), length(windowSize), windowSize(3)]);
window = prod(permute(window, [2, 1, 4, 3]), 4);

end %end function

