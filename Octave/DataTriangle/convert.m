clear
n=[1,3,4,6,7]
for i=1:length(n)
filename=['n',num2str(n(i))]
load(filename)
save('-ascii', ['../DataTriangle3/',filename], 'data')
clear data
end
