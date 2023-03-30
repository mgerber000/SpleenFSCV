function num = monthconvert(str)

switch str
    case 'Jan'
        num = '1';
    case 'Feb'
        num = '2';
    case 'Mar'
        num = '3';
    case 'Apr'
        num = '4';
    case 'May'
        num = '5';
    case 'Jun'
        num = '6';
    case 'Jul'
        num = '7';
    case 'Aug'
        num = '8';
    case 'Sep'
        num = '9';
    case 'Oct'
        num = '10';
    case 'Nov'
        num = '11';
    case 'Dec'
        num = '12';
    otherwise
        num = 'NaN';
end