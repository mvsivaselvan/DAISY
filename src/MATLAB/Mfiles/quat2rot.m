function RR = quat2rot(qq)
normqq = sqrt(qq'*qq);
RR = 1/normqq^2 * ...
    [qq(1)^2+qq(2)^2-qq(3)^2-qq(4)^2 2*(qq(2)*qq(3)-qq(1)*qq(4)) ...
    2*(qq(2)*qq(4)+qq(1)*qq(3)); ...
    2*(qq(2)*qq(3)+qq(1)*qq(4)) qq(1)^2+qq(3)^2-qq(2)^2-qq(4)^2 ...
    2*(qq(3)*qq(4)-qq(1)*qq(2)); ...
    2*(qq(2)*qq(4)-qq(1)*qq(3)) 2*(qq(3)*qq(4)+qq(1)*qq(2)) ...
    qq(1)^2+qq(4)^2-qq(2)^2-qq(3)^2];
