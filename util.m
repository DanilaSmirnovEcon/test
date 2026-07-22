
function u=util(param, a)
if param.nu~=1
    u=(max(a+param.T,param.minValue).^(1-param.nu))./(1-param.nu);
else
    u = log(max(a+param.T, param.minValue));
end
end
