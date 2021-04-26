max_data = ones(1,501);

max_data(120:125) = 0.5;
max_data(127) = 0.5;
max_data(129) = 0.5;
max_data(133:140) = 0.5;

max_data(160:170) = 1;

for time=140:320
    max_data(time) = round(rand())*0.5+0.5;
end
plot(max_data)