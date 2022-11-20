function RMSE= RMSE(yhat, y)
    RMSE = sqrt(mean((y(:) - yhat(:)).^2));  % Root Mean Squared Error

end