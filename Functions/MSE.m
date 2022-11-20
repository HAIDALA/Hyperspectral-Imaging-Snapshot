function MSE= MSE(image_new,lena)
    [M, N] = size(lena);
    error = lena - (image_new);
    MSE = sum(sum(error .* error)) / (M * N);
end