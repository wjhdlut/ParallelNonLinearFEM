# ParallelNonLinearFEM
此项目是本人根据个人兴趣开发的非线性有限元程序。
目前已具备的功能包括：
一、求解流程
1、Newton-Raphson非线性隐式求解，主要处理常规的几何非线性和材料非线性问题。
2、Modified Riks弧长法，主要处理存在Snap through和snap back的后屈曲问题。
3、中心差分法显式求解方法，主要处理高速冲击问题。
二、单元类型
1、三维单元：8节点六面体单元、4节点四面体单元；
2、二维单元：8节点四边形单元、4节点四边形单元、6节点三角形单元、3节点三角形单元、轴对称单元；
3、一维单元：2节点的Timoshenko梁单元、3节点的Timoshenko梁单元、2节点的Kirchhoff梁单元、3节点的Kirchhoff梁单元；
4、单元技术：
   （1）、8节点六面体单元的缩减积分技术，同时支持：Flanagan-Belytschko沙漏控制，主要针对显式动力学问题
三、本构模型
1、各向同性线弹性材料；
2、各向同性弹塑性材料，暂时只支持Von Mises屈服，各向同性硬化模型；
3、各向同行弹性损伤材料，等效应变采用的是Von Mises形式；
四、边界条件
1、约束
   （1）、全约束
   （2）、刚性墙约束
2、载荷
   （1）、指定节点位移
   （2）、指定节点力载荷
   （3）、分布力
五、典型案例
1、Taylor Bar 撞击问题
