% Minimal runtime gate only; no toolbox algorithm or culture model executed.
assert(~usejava('jvm'), 'JVM must remain disabled');
assert(license('test','Symbolic_Toolbox') == 1, 'Symbolic toolbox entitlement unavailable');
disp(version);
syms audit_x positive
assert(isAlways(diff(exp(audit_x),audit_x)==exp(audit_x)));
assert(isAlways(diff(audit_x^audit_x,audit_x)==audit_x^audit_x*(log(audit_x)+1)));
disp('SYMBOLIC_RUNTIME_GATE_COMPLETED');
