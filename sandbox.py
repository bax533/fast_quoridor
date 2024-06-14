import fast_quoridor as qcp
import time

start = time.time()

board = qcp.BoardCpp()

board.CalculatePossibleMoves()
moves_num = board.GetNumOfPossibleMoves()
assert moves_num == 87 # 7*6 * 2 walls + 3 pawn moves
print(moves_num, "number of possible moves")
print(time.time() - start, "s")


print("Moves: ")

moves = board.GetPossibleMoves()
for move in moves:
	move.Print()

print("---")

fence_move = qcp.Move(qcp.M_FENCE, 3, 6, 'h')
pawn_move = qcp.Move(qcp.M_PAWN, 1, 4, 'x')

print("start board:")
board.Print()
print("--")

board.Apply(pawn_move)
print("board after pawn move:")
board.Print()
print("--")

board.Apply(fence_move)
print("board after fence move:")
board.Print()
print("--")
