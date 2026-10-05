"""Build docs/module-c-explainer.pdf: Module C in plain language.

Plain Helvetica only, so no Greek letters or arrows anywhere in the text (the built-in
fonts have no glyphs for them and they come out as black boxes). dG, ddG, ->, x, about.
"""

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase.pdfmetrics import registerFontFamily
from reportlab.platypus import (BaseDocTemplate, Frame, KeepTogether, PageTemplate, Paragraph,
                                Spacer, Table, TableStyle)

OUT = "/Users/vishnut/igem/LOCKR/docs/module-c-explainer.pdf"

INK = colors.HexColor("#1b1b1b")
MUTED = colors.HexColor("#5c5c5c")
ACCENT = colors.HexColor("#1f6f5c")
RULE = colors.HexColor("#d4d4d4")
BAND = colors.HexColor("#f2f5f4")
WARN = colors.HexColor("#8a3d2e")

# Without this, an inline <b> in a Paragraph silently renders at normal weight.
registerFontFamily("Helvetica", normal="Helvetica", bold="Helvetica-Bold",
                   italic="Helvetica-Oblique", boldItalic="Helvetica-BoldOblique")

styles = getSampleStyleSheet()
S = {
    "title": ParagraphStyle("title", parent=styles["Title"], fontName="Helvetica-Bold",
                            fontSize=23, leading=27, textColor=INK, alignment=TA_LEFT,
                            spaceAfter=4),
    "subtitle": ParagraphStyle("subtitle", fontName="Helvetica", fontSize=11.5, leading=15,
                               textColor=MUTED, spaceAfter=16),
    "h1": ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=14, leading=17,
                         textColor=ACCENT, spaceBefore=17, spaceAfter=6),
    "h2": ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=10.8, leading=14,
                         textColor=INK, spaceBefore=10, spaceAfter=3),
    "body": ParagraphStyle("body", fontName="Helvetica", fontSize=10, leading=14.4,
                           textColor=INK, spaceAfter=7),
    "bullet": ParagraphStyle("bullet", fontName="Helvetica", fontSize=10, leading=14.2,
                             textColor=INK, leftIndent=14, bulletIndent=3, spaceAfter=3.5),
    "cell": ParagraphStyle("cell", fontName="Helvetica", fontSize=8.9, leading=11.8,
                           textColor=INK),
    "cellhead": ParagraphStyle("cellhead", fontName="Helvetica-Bold", fontSize=8.9, leading=11.8,
                               textColor=colors.white),
    "mono": ParagraphStyle("mono", fontName="Courier", fontSize=8.6, leading=11.6,
                           textColor=INK, leftIndent=9, spaceBefore=3, spaceAfter=8),
    "quote": ParagraphStyle("quote", fontName="Helvetica-Oblique", fontSize=9.3, leading=12.6,
                            textColor=WARN, leftIndent=12, rightIndent=8, spaceBefore=4,
                            spaceAfter=9),
    "caption": ParagraphStyle("caption", fontName="Helvetica-Oblique", fontSize=8.6, leading=11,
                              textColor=MUTED, spaceBefore=-3, spaceAfter=10),
}


def para(text, style="body"):
    return Paragraph(text, S[style])


def bullets(items):
    return [Paragraph(t, S["bullet"], bulletText="•") for t in items]


def rule(space_before=2, space_after=8):
    t = Table([[""]], colWidths=[6.7 * inch], rowHeights=[0.5])
    t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, -1), 0.6, RULE)]))
    return [Spacer(1, space_before), t, Spacer(1, space_after)]


def table(rows, widths, head=True):
    data = []
    for i, row in enumerate(rows):
        style = "cellhead" if (head and i == 0) else "cell"
        data.append([Paragraph(c, S[style]) for c in row])
    t = Table(data, colWidths=[w * inch for w in widths], hAlign="LEFT", repeatRows=1 if head else 0)
    commands = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
        ("BOX", (0, 0), (-1, -1), 0.5, RULE),
    ]
    if head:
        commands += [("BACKGROUND", (0, 0), (-1, 0), ACCENT),
                     ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, BAND])]
    else:
        commands += [("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.white, BAND])]
    t.setStyle(TableStyle(commands))
    return [t, Spacer(1, 10)]


def decorate(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(0.95 * inch, 0.62 * inch, "LATTE / LOCKR switch design tool")
    canvas.drawRightString(7.65 * inch, 0.62 * inch, "page %d" % doc.page)
    canvas.setStrokeColor(RULE)
    canvas.setLineWidth(0.5)
    canvas.line(0.95 * inch, 0.8 * inch, 7.65 * inch, 0.8 * inch)
    canvas.restoreState()


story = []
A = story.append
E = story.extend

# ---------------------------------------------------------------- title
A(para("Module C: temperature, and working backwards", "title"))
A(para("What the newest part of LATTE does, why it refuses to answer some questions, "
       "and how to read what it gives back.", "subtitle"))
E(rule(0, 4))
A(para("A note on symbols: textbooks write free energy with a capital Greek delta. This guide "
       "writes it dG, and a change in free energy is ddG, which is what the code calls them too. "
       "Energies are in kcal/mol, concentrations in nanomolar (nM), temperatures in kelvin (K). "
       "Room temperature, 25 degrees C, is 298.15 K. Body temperature, 37 degrees C, is 310.15 K.",
       "caption"))

# ---------------------------------------------------------------- short version
A(para("1. The short version", "h1"))
A(para("LATTE already had a calculator that ran one way. You gave it numbers describing a "
       "biosensor and it told you how much brighter that sensor gets when it finds what it is "
       "looking for. That brightness ratio is called <b>fold change</b>. A fold change of 8 means "
       "the sample glows 8 times brighter with the target present than without it."))
A(para("Module C adds two things on top."))
E(bullets([
    "<b>Temperature.</b> The old calculator assumed a sensor behaves the same in a 4 degree C "
    "fridge as in a 37 degree C incubator. It does not. Module C redoes the chemistry at whatever "
    "temperature you actually run the experiment at, and when it does not have enough information "
    "to do that honestly, it says so rather than guessing.",
    "<b>Working backwards.</b> Instead of asking \"what fold change will this design give me\", "
    "you can now ask \"I need a fold change of 8, which designs get me there\". That is the "
    "harder question, and it is the one a lab actually has.",
]))

# ---------------------------------------------------------------- the switch
A(para("2. The sensor, in one paragraph", "h1"))
A(para("The sensor is a protein that works like a lock. A <b>cage</b> sits shut most of the time, "
       "held that way by an arm called the <b>latch</b>. A separate piece called the <b>key</b> "
       "makes light, but only once the cage is open. When the <b>target</b> (a virus protein, say) "
       "binds, it pulls the latch off, the cage opens, the key gets in, and the sample lights up. "
       "Fold change is just the light with target divided by the light without it."))
A(para("How tightly each piece holds on is written as an energy. More negative means a tighter "
       "grip. These are the five numbers Module C can solve for."))
E(table([
    ["Parameter", "What it means in words", "Typical value"],
    ["dG_latch_cage", "How hard the latch grips the cage. Very negative means a cage that stays "
                      "shut, which gives a dark background but is harder to open.", "-4 kcal/mol"],
    ["dG_key_cage", "How hard the light-making key grips the open cage.", "-11 kcal/mol"],
    ["dG_target_binder", "How hard the target grips the part built to catch it.", "-12 kcal/mol"],
    ["key_conc", "How much key you pipette in.", "500 nM"],
    ["target_conc", "How much target is in the sample. This is usually what you are trying to "
                    "measure, not something you choose.", "0.1 to 100 nM"],
], [1.35, 4.05, 1.3]))
A(para("One more number matters and cannot be solved for: <b>pull</b>, how strongly target binding "
       "yanks on the latch. It is set by the shape of the cage and latch, so changing it means "
       "redesigning the protein, not changing a concentration."))

# ---------------------------------------------------------------- part 1
A(para("3. Part 1: temperature", "h1"))
A(para("An energy is not really one number. It is made of a heat part (dH, enthalpy) and a "
       "disorder part (dS, entropy), and those two trade off differently as things warm up. The "
       "equation that describes the trade is the Gibbs-Helmholtz equation, which Module C "
       "implements in full:"))
A(Paragraph("dG(T) = dH_ref - T*dS_ref + dCp*[(T - T_ref) - T*ln(T/T_ref)]", S["mono"]))
A(para("Read it left to right: start from what was measured at a reference temperature, correct "
       "for the new temperature, then apply a second correction using dCp (heat capacity). That "
       "second term is small a few degrees away from the reference and gets important further out, "
       "because proteins change how much heat they soak up as they fold and unfold."))
A(para("Here is the catch, and it is the whole reason this part of the module exists. You can only "
       "run that equation if somebody measured dH. Plenty of published values are just a dG at one "
       "temperature, and from one number there is no honest way to extrapolate. So Module C sorts "
       "whatever you give it into one of three modes and stamps every answer with the mode it used, "
       "so nobody has to guess how much to trust it."))
E(table([
    ["Mode", "What you supplied", "What it does", "How far to trust it"],
    ["<b>full</b>", "dH and dCp", "The whole equation above.", "Good across a wide range."],
    ["<b>vant_hoff</b>", "dH only", "Same equation with dCp taken as zero.",
     "Fine within about 15 K of the reference. Warns past that."],
    ["<b>isothermal</b>", "dG only", "Nothing. Returns the reference value unchanged.",
     "It is not a value at your temperature at all, and it says so."],
], [0.95, 1.15, 2.1, 2.5]))
A(para("That bottom row is the point. Ask for 37 degrees C with only a room-temperature number and "
       "you get the number you started with, plus this:"))
A(para("\"TEMPERATURE SCALING UNAVAILABLE for dG_latch_cage: only a reference dG at 298.15 K was "
       "supplied, so the value at 310.00 K is the UNSCALED 298.15 K value. Supply dH_ref (and "
       "ideally dCp) to extrapolate. Do not read the fold change as a prediction at 310.00 K.\"",
       "quote"))
A(para("It would have been easy to quietly return something that looks like an answer. A number "
       "that looks fine and is wrong is worse than no number, because it ends up in a poster."))
A(para("Two smaller pieces come with this:"))
E(bullets([
    "<b>A dCp estimator.</b> If you know roughly how much surface area gets buried when the two "
    "pieces stick together, there is a helper that turns that into a dCp. It is an estimate from a "
    "published rule of thumb, not a measurement, and any answer that used it carries a line saying "
    "exactly that.",
    "<b>A melting warning.</b> Give it a melting temperature (Tm) and it warns when your working "
    "temperature comes within 10 K of it. Past that point the cage is falling apart on its own, so "
    "the whole open-or-shut picture stops describing reality and the fold change is meaningless "
    "rather than merely imprecise.",
]))
A(para("All three of the older fold-change endpoints now take an optional temperature_K, and it "
       "defaults to 298.15 K. At the default they return exactly what they returned before, to the "
       "last digit. That is checked by a test that compares the two responses field by field."))

# ---------------------------------------------------------------- part 2
A(para("4. Part 2: working backwards", "h1"))
A(para("How much can honestly be said depends entirely on how many parameters you let float. So "
       "there are three different answers rather than one, and the module picks based on the length "
       "of your solve_for list."))
A(para("One free parameter", "h2"))
A(para("There is a single answer and a well-known method (Brent's) finds it quickly. Before "
       "trusting that, the solver checks whether fold change moves in one consistent direction "
       "across your search range. If it ever doubles back, reporting one answer would be hiding "
       "the others, so it reports every solution it can bracket on a fine scan."))
A(para("Two free parameters", "h2"))
A(para("Now there is no single answer. There is a whole curve of combinations that all hit your "
       "target, and the curve is the useful output, so you get it as a list of coordinates ready to "
       "plot. You also get one recommended point picked off that curve: the one sitting furthest "
       "from the edge of the acceptable region."))
A(para("The reason is practical. Think of the acceptable designs as a road. Two designs can both be "
       "on the road, but one is in the middle of the lane and one is half off the shoulder. Pipetting "
       "is never exact and neither is a measured energy, so the middle-of-the-road design is the one "
       "that still works tomorrow. That is what \"most robust\" means here, and it is measured as "
       "distance to the edge in each parameter's own natural units."))
A(para("Three or more free parameters", "h2"))
A(para("The region becomes a shape in three or more dimensions. It cannot be drawn, and mapping it "
       "properly costs far more than it is worth, so the module does not pretend to. Instead you "
       "get the recommended design, and for each parameter separately, how far that one can drift "
       "on its own before the design misses the spec. That per-parameter range is the number you "
       "can actually hand to someone at a bench."))
A(para("When the target is impossible", "h2"))
A(para("This case matters more than it sounds. The chemistry has a hard ceiling: at a pull strength "
       "of 10, fold change cannot get much past about 11, no matter what you do to the energies or "
       "the concentrations. Ask for 1000 and the module does not hand back its closest guess and "
       "hope you do not notice. It reports that the target is out of reach, the best it could "
       "actually do, where that happens, and which parameter is the thing blocking it:"))
# The block and the line explaining it are useless apart, so they move pages together.
A(KeepTogether([
    Paragraph("feasible: false<br/>"
              "message: Target fold change 1000 is not reachable for any dG_latch_cage in [-10, -1].<br/>"
              "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Closest is 11.0 at dG_latch_cage = -10.<br/>"
              "binding_constraint: dG_latch_cage at its lower bound (-10)", S["mono"]),
    para("Read that as: stop tuning the latch, it is already as tight as you allowed, and it is not "
         "the thing standing in your way. Raise the pull, which means redesigning the cage and latch "
         "geometry.", "caption"),
]))

# ---------------------------------------------------------------- outputs
A(para("5. What comes back with every answer", "h1"))
E(table([
    ["Field", "What it tells you"],
    ["values", "The parameter values that hit your target, in kcal/mol and nM."],
    ["fold_change", "What those values actually achieve. Compare it to what you asked for."],
    ["sensitivity", "How much the fold change moves when a parameter moves a little. A big number "
                    "means a touchy design: a small mistake in that parameter shows up as a big "
                    "change in brightness."],
    ["tolerance", "The range each parameter can wander over on its own while the design still "
                  "meets spec. Narrow means fussy."],
    ["warnings", "Everything the temperature layer wanted to tell you, including the refusal to "
                 "extrapolate and the melting-point warning."],
    ["thermo_mode", "full, vant_hoff or isothermal. If several energies are involved, the least "
                    "trustworthy one is reported, because that is the one that limits the answer."],
    ["feasible", "false means the target was out of reach and you should read closest instead of "
                 "looking for a solution."],
], [1.15, 5.55]))

# ---------------------------------------------------------------- part 3
A(para("6. Part 3: turning an energy into DNA you can order", "h1"))
A(para("The solver answers in kcal/mol. It will tell you that you need dG_latch_cage to be -5.5, "
       "which nobody can order. A lab orders DNA coding for specific amino acids, so something has "
       "to bridge the gap."))
A(para("You hand this part your current latch sequence, its current energy, and which positions sit "
       "at the interface. It works out the gap it needs to close (the required ddG), then tries "
       "swapping residues at those positions and scores each swap using the same sequence rules "
       "Module A already uses: acidic residues (D and E) at an interface cost you binding, and "
       "residues that form good helices help hold the latch together. Single swaps and pairs both "
       "get tried, and you get the closest ones back, written the way a lab writes them. D4N means "
       "the D at position 4 becomes an N."))
A(para("Each candidate carries a confidence flag, and the flags are deliberately modest:"))
E(bullets([
    "<b>moderate</b>: the charge rule is doing most of the work, and that rule is calibrated "
    "against a real sequence pair from this project (PfLB-1 and its optimized version).",
    "<b>low</b>: the helix part is doing most of the work, which is a scaling guess rather than a "
    "measurement. Every double mutant is also low, because two changes almost never add up as "
    "neatly as the arithmetic pretends.",
]))
A(para("The response states in plain words that these are heuristic rankings, not predictions, and "
       "that they need checking with Rosetta or FoldX or a plate reader before anyone orders a big "
       "panel of DNA. Ranking is deliberately coarse, to 0.1 kcal/mol, because the estimate is not "
       "better than that. Anything inside one bin counts as a tie, and a single mutation beats a "
       "double in a tie, since one change is cheaper to order and likelier to behave."))

# ---------------------------------------------------------------- tests
A(para("7. How it was checked", "h1"))
A(para("88 new tests, and the suite now stands at 282 passing. Three tests fail, and the same three "
       "were already failing before any of this was written, in unrelated sequence-checking code."))
A(para("The central test is a round trip. Pick random but valid parameters, run the old calculator "
       "forward to get a fold change, hand that fold change to the new solver, and check that the "
       "design it comes back with produces the same fold change. If the two halves ever disagree, "
       "one of them is wrong, and the test says which parameter set broke it."))
A(para("The temperature equation is checked against values worked out by hand at four "
       "temperatures a lab actually uses, using dH = -12.0 kcal/mol, dS = -0.025 kcal/mol/K and "
       "dCp = -0.5 kcal/mol/K:"))
E(table([
    ["Temperature", "Why that temperature", "dG(T), hand-computed"],
    ["277.00 K", "4 degrees C, a cold room or a fridge", "-4.6907 kcal/mol"],
    ["298.15 K", "25 degrees C, the reference point", "-4.5463 kcal/mol"],
    ["310.00 K", "about 37 degrees C, body temperature", "-4.1338 kcal/mol"],
    ["328.00 K", "55 degrees C, hot enough to start worrying about the cage", "-3.0766 kcal/mol"],
], [1.15, 3.55, 2.0]))
A(para("The rest of the suite covers the things most likely to go quietly wrong: unreachable targets "
       "take the reporting path instead of raising an error or inventing a number, isothermal mode "
       "refuses to extrapolate, non-monotonic cases return every root, and 298.15 K reproduces the "
       "old outputs exactly."))

# ---------------------------------------------------------------- limits
A(para("8. What it deliberately does not do", "h1"))
E(bullets([
    "No molecular dynamics, no structure prediction, no GPU. Everything runs on a laptop CPU, and "
    "a full inverse design finishes in well under a tenth of a second.",
    "No screen for it yet. Module C is API only, at POST /api/inverse_design and POST /api/mutations. "
    "The numbers coming back are already in the units the interface uses, so a tab can be added "
    "without touching the maths.",
    "It does not invent thermodynamic data. If you do not have dH, it will tell you what it cannot "
    "do rather than filling the gap with a plausible guess.",
]))

doc = BaseDocTemplate(OUT, pagesize=letter, title="Module C: temperature, and working backwards",
                      author="LATTE", leftMargin=0.95 * inch, rightMargin=0.85 * inch,
                      topMargin=0.85 * inch, bottomMargin=0.95 * inch)
frame = Frame(doc.leftMargin, doc.bottomMargin, 6.7 * inch,
              letter[1] - doc.topMargin - doc.bottomMargin, id="body")
doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=decorate)])
doc.build(story)
print("wrote", OUT)
