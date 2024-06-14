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
	print(move.Printed())

print("---")

fence_move = qcp.Move(qcp.M_FENCE, 3, 6, 'h')
pawn_move = qcp.Move(qcp.M_PAWN, 1, 4, 'x')

print("start board:")
print(board.Printed())
print("--")

board.Apply(pawn_move)
print("board after pawn move:")
print(board.Printed())
print("--")

board.Apply(fence_move)
print("board after fence move:")
print(board.Printed())
print("--")

import fast_quoridor
from fast_quoridor import BoardCpp
from fast_quoridor import Move as MoveCpp

class AgentMC_CPP:
    def __init__(self, n_of_rollouts, rollout_depth):
        self.n_of_rollouts = n_of_rollouts
        self.rollout_depth = rollout_depth

    def best_move(self, board):
        best_score = -self.n_of_rollouts - 1
        best_move = None
        m_it = 0

        board.CalculatePossibleMoves()
        possible_moves = board.GetPossibleMoves()

        for m in possible_moves:
            current_score = 0
            #print(m_it/len(possible_moves))
            for _ in range(self.n_of_rollouts):
                current_score += 1 if board.RolloutMCTS(m, self.rollout_depth) == board.GetTurn() else -1
            if current_score > best_score:
                best_score = current_score
                best_move = m
            m_it += 1

        return best_move

import time

start = time.time()
b = BoardCpp()

MC = AgentMC_CPP(10,5)
print("best move:", MC.best_move(b).Printed())
print("BEST MOVE GENERATION TIME FOR MCTS AGENT 10, 5", time.time()-start)


print("GAME OF TWO AGENTS_MC(10, 5)")
print(b.Printed())

while(b.Winner() == fast_quoridor.INF):
    b.CalculatePossibleMoves()
    m = MC.best_move(b)
    print("WHITE FENCES NUM:", b.GetWhiteFencesNum())
    print("BLACK FENCES NUM:", b.GetBlackFencesNum())
    print("WHITE" if b.GetTurn() == fast_quoridor.WHITE else "BLACK", "MOVE:", m.Printed())
    b.Apply(m)
    
    print(b.Printed())

print("WHITE" if b.Winner() == fast_quoridor.WHITE else "BLACK", "WON THE GAME")



