clear
n=[1, 4, 5, 15]
for i=1:length(n)
filename=['n',num2str(n(i))]
load(filename)
save('-ascii', ['../DataTetrahedral3/',filename], 'data')
save('-mat', ['../DataTetrahedral2/',filename], 'data')
clear data
end
