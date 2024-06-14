#pragma GCC optimize("Ofast,inline,tracer,unroll-loops,vpt,split-loops,unswitch-loops")

// #define DEBUG 1

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <bits/stdc++.h>
#include <math.h>
#include <unordered_set>
#include <unordered_map>
#include <queue>

constexpr int WHITE = 0;
constexpr int BLACK = 1;

constexpr int LEFT = 0;
constexpr int UP = 1;
constexpr int RIGHT = 2;
constexpr int DOWN = 3;

constexpr int ROW_LENGTH = 8;
constexpr int COL_LENGTH = 8;

constexpr uint64_t TOP_ROW = 0xFF00000000000000ull;
constexpr uint64_t BOT_ROW = 0xFFull;
constexpr uint64_t LEFT_COL = 0x8080808080808080ull;
constexpr uint64_t RIGHT_COL = 0x0101010101010101ull;

constexpr uint64_t WHOLE_BOARD = 0xFFFFFFFFFFFFFFFFull;

constexpr int NO_PATH = 99999999;
constexpr int INF = 99999999;




uint64_t SHIFTED(const uint64_t& bitpos, const int& dir)
{
    if(dir == 0)
        return bitpos << 1;
    else if (dir == 1)
        return bitpos << 8;
    else if (dir == 2)
        return bitpos >> 1;
    else if (dir == 3)
        return bitpos >> 8;
    return bitpos;
}

constexpr uint64_t BITPOSITION(const int& row, const int& col)
{
    return 1ull << (row * 8 + col);
}

constexpr int BITINDEX(const uint64_t& onehot) { return __builtin_ctzll(onehot); }

constexpr int BITS_ON(const uint64_t& mask) { return __builtin_popcountll(mask); }

constexpr int GetRow(const uint64_t& pos)
{
    #ifdef DEBUG
        assert(BITS_ON(pos) == 1); 
    #endif
    
        
    return BITINDEX(pos)/COL_LENGTH;
}

constexpr int GetCol(const uint64_t& pos)
{
    #ifdef DEBUG
        assert(BITS_ON(pos) == 1); 
    #endif


    uint64_t COLUMN = RIGHT_COL;
    for(int i=0; i < ROW_LENGTH; i++)
    {
        if(pos & COLUMN) return i;
        COLUMN <<= 1;
    }
    return -1;
}

void print_map(const uint64_t& row=2, const uint64_t& col=2)
{
    uint64_t pos = 1ull << (row * ROW_LENGTH + col);

    for(int y = 7; y >= 0; y--)
    {
        std::cout << y << ": ";
        for(int x = 7; x >= 0; x--)
        {
            uint64_t cur_pos = 1ull << (y * ROW_LENGTH + x);
            std::cout << (pos & cur_pos ? 1 : 0) << " "; 
        }
        std::cout << "\n";
    }
}

void print_bits_index(const uint64_t arg_bits)
{
    uint64_t bits = arg_bits;
    while(bits != 0) // iterate over possible positions
    {
        uint64_t onehot = bits & - bits;
        bits ^= onehot;
        std::cout << BITINDEX(onehot) << " <-bitindex\n";
    }
    std::cout << "-----\n";
}

void print_bits_from_row_col(const int row, const int col)
{
    print_bits_index(BITPOSITION(row, col));
}


constexpr bool IsOutOfBounds(const int& row, const int& col)
{
    return row < 0 || row > 7 || col < 0 || col > 7;
}

bool IsOutOfBoundsWithDir(const int& pos_row, const int& pos_col, const int& dir)
{
    if((pos_col == 7) && (dir == LEFT)) // GOING LEFT ON LEFTMOST COLUMN
        return true;

    if((pos_col == 0) && (dir == RIGHT)) // GOING RIGHT ON RIGHTMOST COLUMN
        return true;
    
    if((pos_row == 0) && (dir == DOWN)) // GOING DOWN ON BOTTOM ROW
        return true;

    if((pos_row == 7) && (dir == UP)) // GOING UP ON TOP ROW
        return true;

    return false; 
}

bool BitsOutOfBoundsWithDir(const uint64_t& bitpos, const int& dir)
{
    #ifdef DEBUG
        assert(BITS_ON(bitpos) == 1); 
    #endif
    
        
    if((bitpos & LEFT_COL) && (dir == LEFT)) // GOING LEFT ON LEFTMOST COLUMN
        return true;

    if((bitpos & RIGHT_COL) && (dir == RIGHT)) // GOING RIGHT ON RIGHTMOST COLUMN
        return true;
    
    if((bitpos & BOT_ROW) && (dir == DOWN))
        return true;

    if((bitpos & TOP_ROW) && (dir == UP))
        return true;

    return false; 
}

constexpr int M_FENCE = 0;
constexpr int M_PAWN = 1;

class Move{
public:
    Move(int t, int r, int c, char o)
        : type(t), row(r), col(c), orientation(o)
        {}

    bool operator==(const Move& other) const
    {
        return type == other.type && row == other.row && col == other.col && orientation == other.orientation;
    }

	std::string Printed() const
	{
		std::string ret = "";
		if(type == M_FENCE) ret += "FENCE ";
		else if(type == M_PAWN) ret += "PAWN ";
		
		ret += std::to_string(row) + " ";
		ret += std::to_string(col) + " ";	
		
		ret += orientation;
		return ret + "\n";
	}

    int type;
    int row;
    int col;
    char orientation;
};

namespace std{
    template <> struct hash<Move>
    {
        size_t operator()(const Move& move) const
        {
            return (hash<int>()(move.type) ^ hash<int>()(move.row) ^ hash<int>()(move.col) ^ hash<char>()(move.orientation));
        }   
    };
}

/*
 * _h_fence_top has bits 1 if ROOT of a horizontal fence TOP of given position in RIGHT direction is on given position
 * _v_fence_right has bits 1 if ROOT of a vertical fence RIGHT of give position in TOP direction is on given position
 *
 */
class BoardCpp {
public:
	BoardCpp()
	: _h_fence_top(0ull),
	  _v_fence_right(0ull),
	  _white_pos((1ull<<4)),
	  _black_pos((1ull<<(7*8+4))),
	  _white_fences(8),
	  _black_fences(8),
	  _who(WHITE),
      _possible_moves(std::unordered_set<Move>())
	{
        srand(time(NULL));
    }

    bool IsReachable(const uint64_t& source_pos, const int& direction,
            const uint64_t& fence_bitpos = 0ull, const char& orientation = 'x') const
    {

        #ifdef DEBUG
            assert(BITS_ON(source_pos) == 1); 
            assert(direction >= 0 && direction <= 3);
            if(orientation != 'x')
            {
                assert(orientation == 'v' || orientation == 'h');
                assert(BITS_ON(fence_bitpos) == 1);
            }
        #endif

        if(BitsOutOfBoundsWithDir(source_pos, direction))
            return false;

        const uint64_t target_pos = SHIFTED(source_pos, direction);

        if(direction == LEFT)
        {
            const uint64_t possible_wall_1 = target_pos;
            const uint64_t possible_wall_2 = !BitsOutOfBoundsWithDir(target_pos, DOWN) ? SHIFTED(target_pos, DOWN) : target_pos;

            if((_v_fence_right & possible_wall_1) | (_v_fence_right & possible_wall_2))
                return false;

            if(orientation == 'v')
            {
                if((possible_wall_1 & fence_bitpos) | (possible_wall_2 & fence_bitpos))
                    return false;
            }
        }
        else if (direction == RIGHT)
        {
            const uint64_t possible_wall_1 = source_pos;
            const uint64_t possible_wall_2 = !BitsOutOfBoundsWithDir(source_pos, DOWN) ? SHIFTED(source_pos, DOWN) : source_pos;

            if((_v_fence_right & possible_wall_1) | (_v_fence_right & possible_wall_2))
                return false;

            if(orientation == 'v')
            {
                if((possible_wall_1 & fence_bitpos) | (possible_wall_2 & fence_bitpos))
                    return false;
            }
        }
        else if (direction == DOWN)
        {
            const uint64_t possible_wall_1 = target_pos;
            const uint64_t possible_wall_2 = !BitsOutOfBoundsWithDir(target_pos, LEFT) ? SHIFTED(target_pos, LEFT) : target_pos;

            if((_h_fence_top & possible_wall_1) | (_h_fence_top & possible_wall_2))
                return false;

            if(orientation == 'h')
            {
                if((possible_wall_1 & fence_bitpos) | (possible_wall_2 & fence_bitpos))
                    return false;
            }
        }
        else if (direction == UP)
        {
            const uint64_t possible_wall_1 = source_pos;
            const uint64_t possible_wall_2 = !BitsOutOfBoundsWithDir(source_pos, LEFT) ? SHIFTED(source_pos, LEFT) : source_pos;

            if((_h_fence_top & possible_wall_1) | (_h_fence_top & possible_wall_2))
                return false;

            if(orientation == 'h')
            {
                if((possible_wall_1 & fence_bitpos) | (possible_wall_2 & fence_bitpos))
                    return false;
            }
        }

        return true;
    }

    int ShortestPathAstar(const uint64_t& pawn_pos, const int& player, const uint64_t& new_fence = 0ull, const char& orientation = 'x') const
    {
        #ifdef DEBUG
            assert(BITS_ON(pawn_pos) == 1); 
            if(orientation != 'x')
            {
                assert(orientation == 'v' || orientation == 'h');
                assert(BITS_ON(new_fence) == 1);
            }
        #endif

        std::priority_queue<std::pair<int, uint64_t> > open;
        std::unordered_map<uint64_t, int> gScore;
        std::unordered_map<uint64_t, int> fScore;
        
        uint64_t bits = WHOLE_BOARD;
        while(bits != 0) // iterate over possible positions
        {
            uint64_t onehot = bits & - bits;
            bits ^= onehot;

            gScore[onehot] = INF;
            fScore[onehot] = INF;
        }

        int end_row = 0;
        if(player == WHITE) end_row = COL_LENGTH-1;
        
        gScore[pawn_pos] = 0;
        fScore[pawn_pos] = abs(end_row - static_cast<int>(BITINDEX(pawn_pos)/COL_LENGTH));
        open.push({0, pawn_pos});        

        const std::unordered_set<int> DIRECTIONS = {LEFT, UP, RIGHT, DOWN};
        while(!open.empty())
        {
            std::pair<int, uint64_t> top = open.top();
            open.pop();
            int current_len = top.first;
            uint64_t current_pos = top.second;

            if((BITINDEX(current_pos) / ROW_LENGTH) == end_row)
                return current_len;

            for(const int& dir : DIRECTIONS)
            {
                if(IsReachable(current_pos, dir, new_fence, orientation))
                {
                    uint64_t neighbour_pos = SHIFTED(current_pos, dir);
                    int new_score = gScore[current_pos] + 1;
                    if(new_score < gScore[neighbour_pos])
                    {
                        gScore[neighbour_pos] = new_score;
                        fScore[neighbour_pos] = new_score + abs(end_row - static_cast<int>(BITINDEX(neighbour_pos)/COL_LENGTH));
                        open.push({fScore[neighbour_pos], neighbour_pos});
                    }
                }
            }
        }

        return NO_PATH;
    }

	bool CanFenceBePlaced(const int& fence_row,const int& fence_col, const char& orientation) const
	{
        #ifdef DEBUG
            // std::cout << fence_row << ", " << fence_col << " CAN FENCE BE PLACED  A\n";
            assert(orientation == 'v' || orientation == 'h');
        #endif
	    if(fence_row >= 7 or fence_row < 0) { return false; }
        if(fence_col > 7 or fence_col <= 0) { return false; }

        uint64_t fence_pos = (1ull << (fence_row * ROW_LENGTH + fence_col));

        if(orientation == 'h')
        {
            if(_v_fence_right & fence_pos)
            { return false; }
            
            if(_h_fence_top & fence_pos || (!IsOutOfBounds(fence_row, fence_col + 1) && (_h_fence_top & (fence_pos << 1))))
            { return false; }            
        }
        else if (orientation == 'v')
        {
            if(_h_fence_top & fence_pos)
            { return false; }
            
            if(_v_fence_right & fence_pos || (!IsOutOfBounds(fence_row - 1, fence_col) && (_v_fence_right & (fence_pos << ROW_LENGTH))))
            { return false; }
        }

        if(ShortestPathAstar(_white_pos, WHITE, fence_pos, orientation) == NO_PATH)
            return false;
        if(ShortestPathAstar(_black_pos, BLACK, fence_pos, orientation) == NO_PATH)
            return false;

        return true;
	
    }

    int WhiteShortestPath() const
    {
        return ShortestPathAstar(_white_pos, WHITE);
    }

    int BlackShortestPath() const
    {
        return ShortestPathAstar(_black_pos, BLACK);
    }

    float Heuristic() const
    {
        float distToEndStraightW = 1.0f;
        float pathToEndW = 1.0f;
        float fencesNumW = 0.75f;

        float whiteValue =
            static_cast<int>(BITINDEX(_white_pos))/8 * distToEndStraightW +
            WhiteShortestPath() * pathToEndW +
            _white_fences * fencesNumW;

        float blackValue =
            (8-static_cast<int>(BITINDEX(_black_pos))/8) * distToEndStraightW +
            BlackShortestPath() * pathToEndW +
            _black_fences * fencesNumW;

        return whiteValue - blackValue;
    }

    int HeuristicWinningPlayer() const
    {
        return Heuristic() > 0 ? WHITE : BLACK;
    }

    int RolloutMCTS(Move move, const int& depth) const
    {
        BoardCpp newBoard = BoardCpp();
        newBoard._h_fence_top = _h_fence_top;
        newBoard._v_fence_right = _v_fence_right;
        newBoard._white_pos = _white_pos;
        newBoard._black_fences = _black_fences;
        newBoard._who = _who;

        int depth_it = 0;
        while(depth_it < depth)
        {
            newBoard.Apply(move);
            newBoard.CalculatePossibleMoves();
            move = newBoard.RandomMove();
            
            depth_it++;
        }

        return newBoard.HeuristicWinningPlayer();
    }

    Move RandomMove() const // assuming possible moves are calculated!
    {
        int idx = std::rand() % _possible_moves.size();
        
        int cnt = 0;
        for(const Move& m : _possible_moves)
        {
            if(cnt == idx)
                return m;
            cnt++;
        }
        return Move(-1, -1, -1, 'x');
    }

    void PlaceFence(const uint64_t& fence_pos, const char orientation, const int player = -99)
    {
        #ifdef DEBUG
            assert(BITS_ON(fence_pos) == 1);
            assert(orientation == 'v' || orientation == 'h');
        #endif
        
            
        if(orientation == 'h')
            _h_fence_top |= fence_pos;
        else if(orientation == 'v')
            _v_fence_right |= fence_pos;

        if(player == WHITE)
            _white_fences--;
        else if (player == BLACK)
            _black_fences--;
    }

    void MovePawn(const uint64_t& pos, int player = -99)
    {
        #ifdef DEBUG
            assert(BITS_ON(pos) == 1);
        #endif
        
        if(player == -99)
            player = _who;

        if(player == WHITE)
            _white_pos = pos;
        else if (player == BLACK)
            _black_pos = pos;
    }

    void Apply(const Move& move)
    {
        if(move.type == M_FENCE)
            PlaceFence(BITPOSITION(move.row, move.col), move.orientation, _who);
        else if (move.type == M_PAWN)
            MovePawn(BITPOSITION(move.row, move.col), _who);

		_who = 1 - _who;
    }

    void CalculatePossibleMoves()
    {
        int num_of_fences = _who == WHITE ? _white_fences : _black_fences;
        uint64_t player_pos = _who == WHITE ? _white_pos : _black_pos;

        _possible_moves.clear();

        if(num_of_fences > 0)
        {
            for(int row = 0; row < COL_LENGTH-1; row++)
            {
                for(int col = 0; col < ROW_LENGTH-1; col++)
                {
                    if(CanFenceBePlaced(row, col, 'h'))
                        _possible_moves.insert(Move(M_FENCE, row, col, 'h'));
                    if(CanFenceBePlaced(row, col, 'v'))
                        _possible_moves.insert(Move(M_FENCE, row, col, 'v'));
                }
            }
        }
        std::unordered_set<int> DIRECTIONS = {LEFT, UP, RIGHT, DOWN};
        for(const int& dir : DIRECTIONS)
        {
            if(!BitsOutOfBoundsWithDir(player_pos, dir) && IsReachable(player_pos, dir)){
                uint64_t new_pos = SHIFTED(player_pos, dir);
                _possible_moves.insert(Move(M_PAWN, GetRow(new_pos), GetCol(new_pos), 'x')); 
            }
        }
    }

    std::unordered_set<Move> GetPossibleMoves() const
    {
        return _possible_moves;
    }

    int GetNumOfPossibleMoves() const
    {
        return _possible_moves.size();
    }

	int GetTurn() const
	{
		return _who;
	}

	int Winner() const
	{
		if(_white_pos & TOP_ROW) return WHITE;
		if(_black_pos & BOT_ROW) return BLACK;
		return INF;
	}

	std::string BoardPrinted() const
	{
		std::string ret = "";

		char white_player_c = 'w';
		char black_player_c = 'b';
		char wall_right_c = '>';
		char wall_top_c = '^';

		for(int y = 7; y >= 0; y--)
		{
			for(int x = 7; x >= 0; x--)
			{
				std::string current_pos_str = "..."; // player v_wall_right, h_wall_top
				const uint64_t cur_pos = BITPOSITION(y, x);
				
				if(cur_pos & _white_pos) current_pos_str[0] = white_player_c;
				if(cur_pos & _black_pos) current_pos_str[0] = black_player_c;
				if(cur_pos & _v_fence_right) current_pos_str[1] = wall_right_c;
				if(cur_pos & _h_fence_top) current_pos_str[2] = wall_top_c;

				ret += current_pos_str + " ";
			}
			ret += "\n";
		}
		return ret;
	}

	uint64_t _h_fence_top, _v_fence_right;
	uint64_t _white_pos, _black_pos;
	int _white_fences, _black_fences;
	int _who;

    std::unordered_set<Move> _possible_moves;
};

PYBIND11_MODULE(fast_quoridor, m){
	m.doc() = "Quoridor board CPP implementation";
    m.attr("LEFT") = pybind11::int_(LEFT);
    m.attr("RIGHT") = pybind11::int_(RIGHT);
    m.attr("UP") = pybind11::int_(UP);
    m.attr("DOWN") = pybind11::int_(DOWN);
    m.attr("WHITE") = pybind11::int_(WHITE);
    m.attr("BLACK") = pybind11::int_(BLACK);
    m.attr("ROW_LENGTH") = pybind11::int_(ROW_LENGTH);
    m.attr("COL_LENGTH") = pybind11::int_(COL_LENGTH);
    m.attr("M_FENCE") = pybind11::int_(M_FENCE);
    m.attr("M_PAWN") = pybind11::int_(M_PAWN);

    pybind11::class_<Move>(m, "Move")
        .def(pybind11::init<int, int, int, char>())
		.def("Printed", &Move::Printed);

	pybind11::class_<BoardCpp>(m, "BoardCpp")
		.def(pybind11::init<>())
		.def("RolloutMCTS", &BoardCpp::RolloutMCTS)
        .def("CalculatePossibleMoves", &BoardCpp::CalculatePossibleMoves)
        .def("RandomMove", &BoardCpp::RandomMove)
        .def("GetPossibleMoves", &BoardCpp::GetPossibleMoves)
        .def("GetNumOfPossibleMoves", &BoardCpp::GetNumOfPossibleMoves)
        .def("Apply", &BoardCpp::Apply)
		.def("GetTurn", &BoardCpp::GetTurn)
		.def("Winner", &BoardCpp::Winner)
        .def("Printed", &BoardCpp::BoardPrinted);
}
