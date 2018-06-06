// Copyright 2018 Haoyu Wang, the University of Chinese Academy of Sciences

import Foundation

// To calculate the minimal polynomial of a square matrix
// @matrix: the matrix you input
// @coefficientMatrix: the coefficient matrix created by μA(A) = 0
// @triangle: "coefficientMatrix" Gauss Eliminated
// @loop, matrixLoop, matrixFinal: variables used in loop
// @count: the number of parameters in the general solution of the homogeneous linear equations
// @rowHasNotZero: count from the last row of triangle, the row which has more than one non-zero element
// @forRank: matrix used in loop to compare the rank with another matrix
// @gridArray: the array which is used to create matrix "forRank"
// @possibleNotZero: the undetermined element which is not linear correlated with undetermined elements whose values are zero
// @zeroCount: the number of undetermined elements whose values are zero
// @gridToSolve: the array which is used to create the augmented matrix which will be used to get the value of each parameter
// @augmentedMatrix: the augmented matrix "gridToSolve" created
// @triangle2: "augmentedMatrix" Gauss Eliminated
// @x: the array which is used to save the value of each parameter
// @sum: variable which is used to calculate the value of each parameter
// @result: the array which is used to save the value of each undetermined element
public func minimalPolynomial(_ matrix: Matrix) -> [Double] {
    
    var triangle: Matrix
    let n = matrix.rows
    var loop = 0
    var matrixLoop = matrix
    var matrixFinal = matrix
    var count = -1
    var rowHasNotZero = 1
    var forRank: Matrix
    var gridArray: [Double] = []
    var possibleNotZero = 0
    var zeroCount = 0
    var gridToSolve: [Double] = []
    let augmentedMatrix: Matrix
    var triangle2: Matrix
    var x: [Double] = []
    var sum: Double = 0
    var result: [Double] = []
    
    assert(matrix.rows == matrix.columns, "Matrix must be square")
    
    var coefficientMatrix = Matrix.init(rows: matrix.rows * matrix.rows, columns: matrix.rows + 1, repeatedValue: 0)
    
    // input values for the last column of coefficientMatrix
    for a in 1...n {
        for j in loop * n + 1...(loop + 1) * n {
            if j - loop * n == a {
                coefficientMatrix[j, n + 1] = 1
            }
            else {
                coefficientMatrix[j, n + 1] = 0
            }
        }
        loop += 1
    }
    
    // input values for other columns of coefficientMatrix
    for i in (1...coefficientMatrix.columns - 1).reversed()  {
        
        if coefficientMatrix.columns - i - 1 > 0 {
            for _ in 1...coefficientMatrix.columns - i - 1 {
                matrixLoop = matrixLoop <*> matrix
                matrixFinal = matrixLoop
            }
        }
        
        matrixLoop = matrix
        
        for j in 1...coefficientMatrix.rows {
            coefficientMatrix[j, i] = matrixFinal.grid[j - 1]
        }
    }
    
    (triangle, _) = coefficientMatrix.gauss()
    
    // processing error
    for i in 0..<triangle.grid.count {
        
        if fabs(triangle.grid[i]) < 0.000000000001 {
            triangle.grid[i] = 0
        }
    }
    
    // @notZeroColumn: the column of triangle which has non-zero element
    // @zeroX: the undetermined element whose value is zero
    // @totalZero: columns the undetermined elements whose values are zero are correponded to
    var notZeroColumn = 0
    var zeroX = 0
    var totalZero: [Int] = []
    
    // get the number of parameters
    for i in (1...triangle.rows).reversed() {
        for j in 1...triangle.columns {
            if triangle[i, j] != 0 {
                count += 1
                notZeroColumn = j
                if totalZero.contains(j) {
                    count -= 1
                }
            }
        }
        if count == 0 {
            count = -1
            if !totalZero.contains(notZeroColumn) {
                zeroX += 1
                totalZero.append(notZeroColumn)
            }
        }
        else if count != -1 {
            rowHasNotZero = i
            break;
        }
    }
    
    // a very rare situation: the number of parameter is 1, and the values of undetermined elements are zero, except the coefficient of the item in minimal polynomial which has the highest degree
    if count == -1 {
        count = 0
    }
    
    // @coefficient2DArray: the two dimensional array which saves the coefficient of each parameter for each undetermined element
    var coefficient2DArray = [[Double]].init(repeating: [0], count: triangle.columns)
    
    // initialize coefficient2DArray
    for i in 1...triangle.columns {
        var coefficientArray:[Double] = []
        if count == 0 {
            count += 1
            for _ in 1...count {
                coefficientArray.append(1)
            }
        }
        else {
            for _ in 1...count {
                coefficientArray.append(0)
            }
        }
        coefficient2DArray[i - 1] = coefficientArray
    }
    
    // continue to initialize coefficient2DArray
    for i in (1...count).reversed() {
        if totalZero.contains(triangle.columns - i - zeroX) {
            coefficient2DArray[triangle.columns - i - zeroX][count - i] = 0
        }
        else {
            coefficient2DArray[triangle.columns - i - zeroX][count - i] = 1
        }
    }
    
    // @j: the column of the first non-zero element in a row
    // @process: variable used in the loop
    var j = triangle.columns - count - zeroX
    var process: Double = 0
    
    // j == 0 being true means that the very rare situation I mentioned before occurred
    // to calculate the coefficient of each parameter for each undetermined element
    if j != 0 {
        for i in (1...rowHasNotZero).reversed() {
            
            for a in 0..<count {
                for b in j...triangle.columns - 1 - zeroX {
                    process += -triangle[i, b + 1] * coefficient2DArray[b][a]
                }
                if triangle[i, j] == 0 {
                    coefficient2DArray[j - 1][a] = 1
                }
                else {
                    coefficient2DArray[j - 1][a] = process / triangle[i, j]
                }
                process = 0
            }
            
            j -= 1
            if j == 0 {
                break
            }
        }
    }
    else {
        for a in 0..<count {
            coefficient2DArray[0][a] = 1
        }
    }
    
    // from row 187 to row 234: to get the array which will be used to create "augmentdeMatrix"
    for a in 0..<count {
        gridArray.append(coefficient2DArray[0][a])
    }
    
    for i in 1..<triangle.columns {
        
        if zeroCount < triangle.columns - i {
            zeroCount = 0
            for j in 1...triangle.columns - i {
                
                for a in 0..<count {
                    gridArray.append(coefficient2DArray[i + j - 1][a])
                }
                
                forRank = Matrix.init(rows: i + 1, columns: count, grid: gridArray)
                
                if forRank.rank() == forRank.rows - 1 { // means linear correlated
                    for _ in 1...count {
                        gridArray.remove(at: gridArray.count - 1)
                    }
                    zeroCount += 1
                    continue
                }
                else if forRank.rank() == forRank.rows { // not linear correlated
                    possibleNotZero = i + j - 1
                    break
                }
            }
        }
        else {
            break
        }
    }
    
    for i in 0...possibleNotZero {
        
        for j in 0..<count {
            gridToSolve.append(coefficient2DArray[i][j])
        }
        
        if i == possibleNotZero {
            gridToSolve.append(1)
        }
        else {
            gridToSolve.append(0)
        }
    }
    
    augmentedMatrix = Matrix.init(rows: possibleNotZero + 1, columns: count + 1, grid: gridToSolve)
    
    (triangle2, _) = augmentedMatrix.gauss()
    
    // processing error
    for i in 0..<triangle2.grid.count {
        
        if fabs(triangle2.grid[i]) < 0.000000000001 {
            triangle2.grid[i] = 0
        }
    }
    
    // from row 248 to row 270: to calculate the value of each parameter
    var i = triangle2.rows - 1
    let n1 = triangle2.rows
    
    for j in 0..<n1 {
        if triangle2[n1 - j, triangle2.columns - 1] == 0 {
            i -= 1
            continue
        }
        else {
            x.append(triangle2[n1 - j, triangle2.columns] / triangle2[n1 - j, triangle2.columns - 1])
            break
        }
    }
    
    while i >= 1 {
        sum = 0
        for j in i + 1...(triangle2.columns - 1) {
            sum += triangle2[i, j] * x[j - i - 1]
        }
        x.insert(((triangle2[i, triangle2.columns] - sum) / triangle2[i, i]), at: 0)
        i -= 1
    }
    
    // @processResult: variable used in the loop
    var processResult: Double = 0
    
    // to calculate the value of each undetermined element x
    for i in 0..<triangle.columns {
        for j in 0..<count {
            processResult += x[j] * coefficient2DArray[i][j]
        }
        result.append(processResult)
        processResult = 0
    }
    
    // processing error
    for i in 0..<result.count {
        if fabs(result[i]) < 0.0000000001 {
            result[i] = 0
        }
    }
    
    // processing error
    // @resultToString: variable used to process error
    var resultToString:String
    
    for i in 0..<result.count {
        resultToString = String(result[i])
        while resultToString.count > 10 {
            resultToString.remove(at: resultToString.index(before: resultToString.endIndex))
        }
        result[i] = Double(resultToString)!
    }
    
    // processing error
    for i in 0..<result.count {
        if fabs(result[i]) - fabs(Double(Int(result[i]))) < 0.000001 || fabs(result[i]) - (fabs(Double(Int(result[i]))) + 1) < 0.000001 || fabs(result[i]) - (fabs(Double(Int(result[i]))) - 1) < 0.000001 {
            if fabs(result[i]) - fabs(Double(Int(result[i]))) < 0.000001 && fabs(result[i]) - fabs(Double(Int(result[i]))) > -0.000001 {
                result[i] = Double(Int(result[i]))
            }
            else if fabs(result[i]) - (fabs(Double(Int(result[i]))) + 1) < 0.000001 && fabs(result[i]) - (fabs(Double(Int(result[i]))) + 1) > -0.000001 {
                if result[i] < 0 {
                    result[i] = Double(Int(result[i])) - 1
                }
                else {
                    result[i] = Double(Int(result[i])) + 1
                }
            }
            else if fabs(result[i]) - (fabs(Double(Int(result[i]))) - 1) < 0.000001 && fabs(result[i]) - (fabs(Double(Int(result[i]))) - 1) > -0.000001 {
                result[i] = Double(Int(result[i])) + 1
            }
        }
    }
    
    return result  // "result" is the array which contains the value of each undetermined element x, or in other words,
                   // the coefficient of each item in the minimal polynomial of "matrix"
}

