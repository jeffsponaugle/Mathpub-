// CPDivide.java - perft with the Chesspresso library, split by first move.
//
// Mirrors Sean A. Irvine's A048987 recursion in joeis (which produced the OEIS terms of
// A285873-A285878): Position.getAllMoves(), doMove(), undoMove(), counting every game that
// reaches the target ply. Each first move's subtree runs on its own thread.
//
// Usage: java -cp build CPDivide "<FEN>" depth threads [move ...]
//   Moves (from and to squares, e.g. e2e4, plus a promotion letter) are played first; the
//   resulting position's FEN is printed, then perft(depth) from it, split by move.
import chesspresso.Chess;
import chesspresso.move.Move;
import chesspresso.position.FEN;
import chesspresso.position.Position;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class CPDivide {
  static String str(short m) {
    String s = Chess.sqiToStr(Move.getFromSqi(m)) + Chess.sqiToStr(Move.getToSqi(m));
    if (Move.isPromotion(m)) s += Character.toLowerCase(Chess.pieceToChar(Move.getPromotionPiece(m)));
    return s;
  }

  static long perft(Position p, int depth) throws Exception {
    if (depth == 0) return 1;
    long n = 0;
    for (short m : p.getAllMoves()) {
      p.doMove(m);
      n += perft(p, depth - 1);
      p.undoMove();
    }
    return n;
  }

  static Position at(String fen, List<String> moves) throws Exception {
    Position p = new Position(fen, true);
    for (String mv : moves) {
      Short found = null;
      for (short m : p.getAllMoves()) if (str(m).equals(mv)) found = m;
      if (found == null) throw new IllegalArgumentException("no move " + mv + " in " + FEN.getFEN(p));
      p.doMove(found);
    }
    return p;
  }

  public static void main(String[] a) throws Exception {
    final String fen = a[0];
    final int depth = Integer.parseInt(a[1]);
    final int threads = Integer.parseInt(a[2]);
    final List<String> pre = Arrays.asList(a).subList(3, a.length);
    final Position root = at(fen, pre);
    System.out.println("fen: " + FEN.getFEN(root));
    final short[] moves = root.getAllMoves();
    final ExecutorService ex = Executors.newFixedThreadPool(threads);
    final List<Future<Long>> fs = new ArrayList<>();
    for (short m : moves) {
      final List<String> line = new ArrayList<>(pre);
      line.add(str(m));
      fs.add(ex.submit(() -> perft(at(fen, line), depth - 1)));
    }
    long total = 0;
    for (int i = 0; i < moves.length; i++) {
      final long c = fs.get(i).get();
      total += c;
      System.out.println(str(moves[i]) + ": " + c);
    }
    ex.shutdown();
    System.out.println("total: " + total);
  }
}
