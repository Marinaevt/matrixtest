�������� dll ���������� ImprovedSystem ������ 1 �� 4 ������� 2006 �.

���������� �������� �������:

������� - ImprovedGaussSystem
��������: ������������� ��� ������� ���� ������� ������
����������: 
1. ����� � �������� ������ ���� Double - 8 ���� 
2. � ���� ������ ��� �� ������������ ������������� ����������
���������: 
function ImprovedGaussSystem(A,C,B,PB:Pointer;NN,Q:Integer):integer;cdecl;
�������:
A - ������ ������� A �� ����� Double(8 ����)(����� ����� ����) ,������ ������ ���� ������ NN*Q*8 ���� � ���������
��������� �����(������� �������� �����(������������ ������)) �.�. ������� �� �������.
C - ������ ������� � �� ����� Double(8 ����)(������ ����� ����) ,������ ������ ���� ������ NN*8 ���� � ��������� ������ ������ ����� ����
B - ������ ������ ���� ������ Q*8 ����, ���������� �� ����� ��������
PB - ������ ������ ���� ������ Q ����, ���������� �� ����� ��������
NN - ����� ����� ������� A
Q - ����� ��������(������ ������(�������� ������ �����)) ������� �
���������: 
������� ���������� ����� ���� ������, � ������ ������ ��� ������ 1
� ������� C ��������� ������� ����
� �������� A,B � PB ��������� ����������� ����� (������� �������� �����)
