import numpy as np



def is_diagonally_dominant(mat):
    if mat is None: #empty matrix
        return False

    d = np.diag(np.abs(mat))  # Find diagonal coefficients
    s = np.sum(np.abs(mat), axis=1) - d  # Find row sum without diagonal
    return np.all(d > s)

def dominantDiagonalTry(mat):
    n =len(mat)
    dom = [-1]*n #working with index so -1 is safe
    matRes = [row.copy() for row in mat]
    
    for i in range (n):
        row = mat[i]
        row_sum = sum(abs(a) for a in row) #sum the row
        for j in range(n):
            if abs(row[j])>(row_sum-abs(row[j])):#dominance for diagonal
                dom[i] = j
                break #found dominance we can stop for row
    if -1 in dom: #cant make dominant diagonal 
        print("no dominant diagonal that we can find")
        return mat
    for i in range(n):#change rows for dominant diagonal
        j = dom[i]
        matRes[j] = mat[i] #change rows
    return matRes
    
def jacobi(A,b,X0arr,ep,N):
    n = len(A)
    k = 1
    

    mass = "iterations\t" + "\t".join([f"x[{pla}]" for pla in range(n)])
    print(mass)
    while k <= N:
        ans= np.zeros(n,dtype = np.double)#array of double with 0 for ansers
        for i in range(n):
            sumRow = 0
            for j in range(n):
                if j != i:#sum all except x in that row
                    sumRow = sumRow + (A[i][j] *X0arr[j]) #sum by x vactor and values
            ans[i] = (b[i]- sumRow)/A[i][i] #jacobian for x number i of vector
            
        print(f"{k}\t" + "\t".join([f"{xi:.8f}" for xi in ans]))
    
        #need to add chaking for exiding for exiting 
        normAns = np.linalg.norm(ans, np.inf)
        normX = np.linalg.norm(X0arr, np.inf)
        if  abs(normAns - normX) < ep:#norm in epsilon of anser
            print("found")
            if(is_diagonally_dominant(A)):
                print("even that we have no dominant diagonal we got anser")
            return tuple(ans)
        k=k+1
        X0arr = ans.copy()#move new values to old 
    print ("anser with epsilon from resulf not found exceed number of iterations")
    return tuple(ans)

def gauss_seidel(A,b,X0arr,ep,N):
    n = len (A)
    k =1
    
    mass ="iterations\t" + "\t".join([f"x[{pla}]" for pla in range(n)])
    print (mass)
    
    while k <=N:
        ans = np.zeros(n,dtype = np.double)
        for i in range(n):
            sumRow = 0
            for j in range(n):
                if j != i:
                    sumRow += A[i][j]*ans[j] #all the diffrense from jacobi
            ans[i] = (b[i] - sumRow)/A[i][i]
        
        print(f"{k}\t" + "\t".join([f"{xi:.8f}" for xi in ans]))
        normAns = np.linalg.norm(ans, np.inf)
        normX = np.linalg.norm(X0arr, np.inf)
        if  abs(normAns - normX) < ep:#norm in epsilon of anser
            print("found")
            if(is_diagonally_dominant(A)):
                print("even that we have no dominant diagonal we got anser")
            return tuple(ans)
        k=k+1
        X0arr = ans.copy()#move new values to old 
    print("Solution did not converge within the maximum number of iterations.")
    return tuple(ans)
        

matrixA = [[4,2,0],[2,10,4],[0,4,5]]
matrixB = [[0.5,8,0],[2,1,4],[0,4,2]] #matrix for chaking
vectorB=[[2],[6],[5]]
erortol=0.00001 #epsilon
N = 200 #number of iterations
x_vec = np.zeros_like(vectorB, dtype=np.double)

'''
#print(is_diagonally_dominant(matrixA))
#
if (is_diagonally_dominant(matrixB)):
    print("true")
else:
    print("False")
    matrixB = dominantDiagonalTry(matrixB)
    print(matrixB)'''
    
print(matrixB)
user_choice = input("Please choose an option 1 for jacobi, 2 for gauss seidel: ")
# Check the user's choice and display the corresponding option
if user_choice == "1":
    print("You chose Jacobi iterations ")
    sol = jacobi(matrixA,vectorB,x_vec,erortol,N)
elif user_choice == "2":
    print("You chose gauss seidel algorithm")
    sol = gauss_seidel(matrixA,vectorB,x_vec,erortol,N)
