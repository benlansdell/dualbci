%regression_nev(nevfile, fn_out, kernellength, binsize, sigma, offset)
nevfile = './testdata/20130117SpankyUtah001.nev';
kl = 6;
verbosity = 1;
offset = 0;

%Sample rate 20Hz, sigma 5;
bs = 0.05;
sigma_fr = 5;
sigma_trq = 10;
fn = ['./worksheets/preprocessing/01172013_sr_' num2str(1/bs) '_sigma_fr' num2str(sigma_fr) '_sigma_trq_' num2str(sigma_trq) '/20130117SpankyUtah001'];
regression_nev(nevfile, fn, kl, bs, sigma_fr, sigma_trq, offset, verbosity);

%Sample rate 100Hz, sigma 5;
bs = 0.01;
sigma_fr = 5;
sigma_trq = 10;
fn = ['./worksheets/preprocessing/01172013_sr_' num2str(1/bs) '_sigma_fr' num2str(sigma_fr) '_sigma_trq_' num2str(sigma_trq) '/20130117SpankyUtah001'];
regression_nev(nevfile, fn, kl, bs, sigma_fr, sigma_trq, offset, verbosity);

%Sample rate 500Hz, sigma 5;
bs = 0.002;
sigma_fr = 5;
sigma_trq = 10;
fn = ['./worksheets/preprocessing/01172013_sr_' num2str(1/bs) '_sigma_fr' num2str(sigma_fr) '_sigma_trq_' num2str(sigma_trq) '/20130117SpankyUtah001'];
regression_nev(nevfile, fn, kl, bs, sigma_fr, sigma_trq, offset, verbosity);

%Sample rate 100Hz, sigma 50;
bs = 0.01;
sigma_fr = 5;
sigma_trq = 50;
fn = ['./worksheets/preprocessing/01172013_sr_' num2str(1/bs) '_sigma_fr' num2str(sigma_fr) '_sigma_trq_' num2str(sigma_trq) '/20130117SpankyUtah001'];
regression_nev(nevfile, fn, kl, bs, sigma_fr, sigma_trq, offset, verbosity);

%Sample rate 100Hz, sigma 20;
bs = 0.01;
sigma_fr = 5;
sigma_trq = 20;
fn = ['./worksheets/preprocessing/01172013_sr_' num2str(1/bs) '_sigma_fr' num2str(sigma_fr) '_sigma_trq_' num2str(sigma_trq) '/20130117SpankyUtah001'];
regression_nev(nevfile, fn, kl, bs, sigma_fr, sigma_trq, offset, verbosity);
