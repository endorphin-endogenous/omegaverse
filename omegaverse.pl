#!/usr/bin/perl
#use warnings;
use Data::Dumper;
use Getopt::Std;
use Config::Tiny;
use Time::HiRes;
#動作要件：perl 5.20.2 および上記モジュール
#引数指定###############################################################
#-i:初期個体数 default:10000
#-a, -b, -o :初期ω遺伝子頻度(整数比で指定可, default: a:b:o = 2:4:1)
#-g:ω遺伝子モデルの指定(default:GeneModel1, サブルーチンOmegatypeを参照)
#-m:繁殖モデルの指定(default:MateModel1, サブルーチンMatingを参照)
#-o:出力ファイルの指定
#-y:経過年数の指定(default:1)

my %opts = ();
getopts("i:a:b:o:g:f:y:m:", \%opts);
my $StartIndividuals = $opts{i} //= 10000;
my $Alphas = $opts{a} //= 2; 
my $Betas = $opts{b} //= 4;
my $Omegas = $opts{o} //= 1;
my $Model = $opts{g} //= "GeneModel1" ;
my $MateModel = $opts{m} //= "MateModel1" ;
my $output = $opts{f};
my $year = $opts{y} //= "1";

#割合に変換
my $all = $Alphas + $Betas + $Omegas;
$Alphas = $Alphas / $all * 100;
$Betas = $Betas / $all * 100;
$Omegas = $Omegas / $all * 100;

#同一ディレクトリ内の設定ファイル(omegaverse.ini)読み込み###################
#1.omegaverse.ini:初期集団の年齢構成および繁殖設定ファイル(引数に組み込み予定)
#2.primes.txt    :素数列(各個体の識別および血縁度判定用の数値)
#3.mortality.txt :死力表(各個体の死亡判定用)
#上２つは現時点では無駄な処理なので今後修正予定

my $cabin = $0;
$cabin =~ s/omegaverse\.pl$//;
my $inicabin = $cabin . "omegaverse.ini";
my $init;
unless(-e $inicabin){
	print "[ERROR] Not Found '$inicabin' in '$cabin'\n";
	exit;
}else{
    $init = Config::Tiny -> read("$inicabin");
}

#ID用の素数列読み込み
my $prime = $0;
$prime =~ s/omegaverse\.pl$//;
$prime = $prime . "primes.txt";
open(PR, "<$prime") or die;
my @primes = ();
my $c = 0;
while(<PR>){
    $c++;
    chomp;
    push(@primes, $_);
    last if $c == $StartIndividuals + 1;
}

#死亡判定用の死力表読み込み(平成２２年度完全生命表)
my $mortality = $0;
$mortality =~ s/omegaverse\.pl$/mortality.txt/;
open(MOR, "<$mortality") or die;
my %death_XX;
my %death_XY;
while(<MOR>){
    s/\r|\n//g;
    my @line = split("\t", $_);
    $death_XY{$line[0]} = $line[1];
    $death_XX{$line[0]} = $line[2];
}
close(MOR);

#初期集団生成#################################################
#各個体のプロファイル作成
my $Omegaverse;
#個体ごとの処理
for($i = 1; $i <= $StartIndividuals; $i++){
    my %Individual;
    $Individual{ID} = $primes[$i];
    $Individual{ParentID} = "0"; #初期親

    #ヒト性別遺伝子型
    my $XYtype = rand(100);
    if ($XYtype > 50){
        $Individual{XYtype} = "XX";
        $Individual{XYcode} = "1";
    }else{
        $Individual{XYtype} = "XY";
        $Individual{XYcode} = "0";
    };

    #オメガ遺伝子型
    my @ABOgene = ();
    for($j = 0; $j < 2; $j++){
        my $gene = rand(100);
        if($gene <= $Alphas){
            push(@ABOgene, "A"); #alpha
        }elsif($gene <= $Betas + $Alphas){
            push(@ABOgene, "B"); #beta
        }else{
            push(@ABOgene, "O"); #omega
        };
    }
    @ABOgene = sort{$a cmp $b} @ABOgene;
    $Individual{ABO} = join("", @ABOgene);

    #オメガ表現型(サブルーチン'Omegatype'で遺伝子モデルごとに決定)
    my %code = ("Alpha" => 2, "Beta" => 1, "Omega" => 0);
    my $omegatype = &Omegatype(\%Individual);
    $Individual{omegatype} = $omegatype;
    $Individual{omegacode} = $code{$omegatype};

    #年齢(世界統計2000年時点を採用、設定ファイルに記載)
    my $age = rand(100);
    if($age < $init->{population}->{under15}){ #15歳未満
        $age = int(rand(15)+0.5);
    }elsif($age >= (100-$init->{population}->{over65})){ #65歳以上
        my $lifespan = $init->{population}->{lifespan} - 65;
        $age = int(rand($lifespan)+0.5) + 65;
    }else{ #15~65歳
        my $lifespan = 50;
        $age = int(rand($lifespan)+0.5) + 15;
    }
    $Individual{age} = $age;

    #つがい(初期は無しとする？)
    my $matepair = 0;
    $Individual{matepair} = $matepair;
    $Individual{matenum} = $matepair;

    #状態(交配可能か、そうでないか)
    my $status = &Matestatus(\%Individual); #交配状態判定用ルーチンへ送る
    $Individual{status} = $status;

    #個体群に代入
    $Omegaverse->{$i} =  \%Individual;
}

#初期状態出力
print "---Inital State---\n";
&IndvidualStat();
&Outstat("0","stat");
print "------------------\n";
#print Dumper $Omegaverse;

#メイン処理##########################################################
#-yでの指定年数間、多次元配列Omegaverseに対してmating-birth-deathの処理を反復
#IndividualStatで画面に統計出力、Outstatでtsvファイルに出力
#
#1.Mating:無作為に2個体抽出→交配状態・遺伝子型・出生率調整に基づきつがい判定
#2.Birth :つがい形成した個体の情報から子のデータ生成、Omegaverseに追加
#3.Death :年齢に基づき死亡判定、生存個体のみをOmegaverseに追加、
#         Beta個体の死亡したつがい情報を削除(alhpa, omegaは固定)
#
#1.2.は人口千人比の婚姻率が5パーセント・出生率7パーセント前後になる程度に調整(omegaverse.ini参照)
for($k = 1; $k <= $year; $k++){
    #結婚処理
    $Omegaverse = &Mating($Omegaverse);
    #繁殖処理
    &Birth($Omegaverse);
    #死亡・生存処理
    $Omegaverse = &Death($Omegaverse);
    my $chr = length($k);
    $chr = "0" x (4 - $chr) . $k;
    print "---year:$chr------\n";
    &IndvidualStat();
    #経過出力
    &Outstat($k,"stat");
}

#メイン処理終了, 以下サブルーチン#####################################

#個体群の統計値モニタリング
sub IndvidualStat {
    #カウント対象のkeyを作成(デフォルトはXYとABO)
    my $countkey = $_[0] //= "ABO|XYtype|omegatype";
    
    #各遺伝子型のカウント
    my %stats;
    foreach my $id (keys($Omegaverse)){
        foreach my $key (keys(%$Omegaverse{$id})){
           next unless ($key =~ /$countkey/);
           my $type = $$Omegaverse{$id}{$key};
           $stats{$type} //= 0;
           my $num = $stats{$type} + 1;
           $stats{$type} = $num;
        }
    }

    if($countkey =~ /ABO\|XYtype/){
        #モデルごとに表現型比、遺伝子型比を計算
        my $total = $stats{XX} + $stats{XY};
        print "Population size: $total\n";
        print "Sex:             ";
        #出力
        #sex
        foreach my $keys (sort keys(%stats)){
            next if ($keys =~ /A|B|O/);
            my $per = $stats{$keys} / $total * 100;
            $per = sprintf("%.2f", $per); #小数点下2ケタまで
            print "$keys: $per%\t";
        }
        print "\n";

        #omegatype
        print "Omegatype:       ";
        foreach my $keys (sort keys(%stats)){
            next if ($keys =~ /^..$/);
            my $per = $stats{$keys} / $total * 100;
            $per = sprintf("%.2f", $per); #小数点下2ケタまで
            print "$keys: $per%\t";
        }
        print "\n";

        #genotype
        print "Genotype:        ";
        foreach my $keys (sort keys(%stats)){
            next if ($keys =~ /XX|YY|XY|.{3,5}/);
            my $gene = $stats{$keys} / $total * 100;
            $gene = sprintf("%.2f", $gene);
            $stats{$keys} = $gene;
            print "$keys: $stats{$keys}% ";
            }
        print "\n";
    }
    #年齢分布
    my %agehist;
    if($countkey =~ /age/){
        my $total = 0;
    }
    #婚姻率分布
    if($countkey =~ /status/){
        print Dumper %stats;
    }
    if($countkey =~ /matepair/){
        print Dumper %stats;
    }

}
#ファイル出力
sub Outstat {
        my $mode = $_[1] //= "stat";
        my $year = $_[0] //= "1";
        my $style = $_[2] //= "per";

        #全個体のデータ出力        
        if($mode eq "mass"){
            open(OUT, ">>$output");
            #header
            my @head = sort (keys(%$Omegaverse{1}));
            print OUT "#";
            foreach my $head (@head){
             print OUT "$head\t";
            }
            print OUT "year\n";
            foreach my $id (sort keys(%$Omegaverse)){
                foreach my $elem (sort keys(%$Omegaverse{$id})){
                    print OUT "$$Omegaverse{$id}{$elem}\t";
                }
                print OUT "$year\n";
            }

        #現時点の個体群の統計出力
        }elsif($mode eq "stat"){
            open(OUT, "+>>$output");
            #集計処理
            my $stats;
            foreach my $id (keys($Omegaverse)){
                foreach my $key (keys(%$Omegaverse{$id})){
                next if $key =~ /ID|age|mate|code/;
                my $type = $$Omegaverse{$id}{$key};
                $$stats{$key}{$type} //= 0;
                my $num = $$stats{$key}{$type} + 1;
                $$stats{$key}{$type} = $num;
                }
            }
            #ヘッダ生成
            #print Dumper $stats;
            if ($year eq 0){
                #引数の記述
                print OUT "#./omegaversepl -i $StartIndividuals -a $Alphas -b $Betas -c $Omegas -g $Model -m $MateModel -y $opts{y}\n";
                #ヘッダ
                print OUT "#year\tTotal\t";
                foreach my $key (sort(keys($stats))){
                    foreach my $stat (sort(keys(%$stats{$key}))){
                        print OUT "$stat\t";
                    }
                }
                print OUT "\n";
            }
            #出力(デフォルトはパーセンテージ)
            my @population = keys($Omegaverse);
            my $total = @population;
            print OUT "$year\t$total\t";
            foreach my $key (sort(keys($stats))){
                foreach my $stat (sort(keys(%$stats{$key}))){
                    my $out = $$stats{$key}{$stat};
                    if ($style eq "per"){
                        my $per = $out / $total * 100;
                        $out = sprintf("%.2f", $per);
                    }
                    print OUT "$out\t";
                }
            }
            print OUT "\n";
        }
    #close(OUT);
}
#オメガ表現型判定(既存・新規モデルごとに編集)
sub Omegatype {
    my $in = $_[0];
    my %input = %$in;
    my $type;
    #GeneModel1 (A>B>O, XY独立)
    if($Model eq "GeneModel1"){
        if($input{ABO} =~ /A/){
            $type = "Alpha";
        }elsif($input{ABO} =~ /B/){
            $type = "Beta";
        }else{
            $type = "Omega";
        }
    #GeneModel2 (A＝O > B, AO = B)
    }elsif($Model eq "GeneModel2"){
        if($input{ABO} =~ /A/){
            $type = "Alpha";
            $type = "Beta" if $input{ABO} eq "AO";
        }elsif($input{ABO} =~ /O/){
            $type = "Omega";
        }else{
            $type = "Beta";
        }
    }
    return($type);
}
#繁殖状態判定
sub Matestatus {
    my $in = $_[0];
    my %input = %$in;
    my $status = 1;
    #初潮前は繁殖不可(仮：男女共通にしておく)
    $status = 0 if ($input{age} < $init->{mating}->{menophania});
    #閉経後は繁殖不可
    $status = 0 if ($input{omegatype} =~ /O/ and $input{age} > $init->{mating}->{menopause});
    $status = 0 if ($input{XYtype} =~ /XX/ and $input{age} > $init->{mating}->{menopause});
    $status = 0 if ($input{XYtype} =~ /XY/ and $input{age} > $init->{mating}->{menopause});

    return($status);
}

#繁殖プロセス(繁殖可不可の判定→子の生成)
sub Mating {
    my $in = $_[0];
    my @pop = sort keys(%$in);
    my $population = @pop;
    my $mate_num = $population * $init->{mating}->{mating_per};
    my $mate_success = 0;
    my $birth_sucsess = 0;
    my $test_num = 0;
    #mate pair生成
    for($j = 1; $j <= int($mate_num); $j++){
        my $mate1 = int(rand($population)+0.5);
        my $mate2 = int(rand($population)+0.5);
        #自殖(同一ID)の場合は再代入
        while($mate1 == $mate2){
            $mate2 = int(rand($population)+0.5);
        }

        #statusに基づき交配可能判定
        next if $$in{$mate1}{status} eq 0;
        next if $$in{$mate2}{status} eq 0;

        #性別に基づいた子の生成可能判定
        my $materesult_omega = $$in{$mate1}{omegacode} + $$in{$mate2}{omegacode};
        my $materesult_XY = $$in{$mate1}{XYcode} + $$in{$mate2}{XYcode};
       
        #$test_num++; #debug
        #MateModel1(XX-XY, XX-XXかつ一方がA, XY-XYかつ一方がO, 同性でA-A, O-Oは除外)
        if($MateModel eq "MateModel1"){
            #同性かつアルファまたはオメガが無いケース,両方がアルファまたはオメガのケースを除外
            if($materesult_XY eq 0){ #XY-XYの場合
                #同性betaを除外
                next if($$in{$mate1}{omegacode} * $$in{$mate2}{omegacode} eq 1);
                #一方のみがOでないと子供は生まれない
                next if($materesult_omega =~ /0|3|4/);
            }elsif($materesult_XY eq 2){  #XX-XXの場合
                #同性betaを除外
                next if($$in{$mate1}{omegacode} * $$in{$mate2}{omegacode} eq 1);
                #一方がAでないと子供は生まれない
                next if ($materesult_omega =~ /0|1|4/);
            }
        #MateModel2(男性O-O, 女性A-Aを許容)
        }elsif($MateModel eq "MateModel2"){
            #同性かつアルファまたはオメガが無いケース,両方がアルファまたはオメガのケースを除外
            if($materesult_XY eq 0){ #XY-XYの場合
                #同性betaを除外
                next if($$in{$mate1}{omegacode} * $$in{$mate2}{omegacode} eq 1);
                #少なくとも一方がOでなければ子供は生まれない
                next if($materesult_omega =~ /3|4/);
            }elsif($materesult_XY eq 2){  #XX-XXの場合
                #同性betaを除外
                next if($$in{$mate1}{omegacode} * $$in{$mate2}{omegacode} eq 1);
                #少なくとも一方がAでないと子供は生まれない
                next if ($materesult_omega =~ /0|1/);
            }
        #エラー処理
        }else{
            print "[ERROR]:'$MateModel' is not supported mating model\n";
            exit;
        }
        $mate_success++;

        #つがいの判定(A-O, およびB-B間ではつがい形成)
        if($$in{$mate1}{omegacode} + $$in{$mate1}{omegacode} eq 2){
            $$in{$mate1}{matepair} = $$in{$mate2}{ID};
            $$in{$mate1}{matenum} = $mate2;
            $$in{$mate2}{matepair} = $$in{$mate1}{ID};
            $$in{$mate2}{matenum} = $mate1;
        }
        #血縁度判定:素因数分解して共通する親の個数が2個
        #17/9/6 血縁度処理はいったん放置：
        #my $kindred = &Kindred($$in{mate1}{ID}, $$in{mate1}{ID}, $population); #サブルーチンへ投げる
        #next if $kindred < 3;
        
        #出生率調整
        my $birth_fate = rand(1);
        next if $birth_fate > $init->{mating}->{birth_rate};
        #$birth_success++;
    }
    #結婚率出力
    my $mating_rate = $mate_success / $population * 1000;
    $mating_rate = sprintf("%.2f", $mating_rate);
    print "Mating rate: $mating_rate\n";
    return($in);

}

#血縁度処理（与えられた2数を素因数分解して、共通する因数の個数を返す）
sub Kindred {
#開発中止
}

#つがいの情報から子供の情報生成、omegaverseに追加
sub Birth {
    my $in = $_[0];
    my @pop = sort keys(%$in);
    my $population = @pop;
    foreach my $id (sort keys(%$in)){
        next if $$in{$id}{matepair} eq 0; #非婚は除外
        next if $$in{$id}{status} eq 2;   #出産済の場合は除外
        #親の情報取り込み
        my $parentid_1 = $$in{$id}{ID};
        my $pair_id = $$in{$id}{matenum};
        my $parentid_2 = $$in{$pair_id}{ID};
        #親の情報変更
        $$in{$id}{status} = 2;
        $$in{$pair_id}{status} = 2;

        #子供の情報生成 parentIDは血縁度の繁殖制限使わない限り無駄そう
        my %child;
        $child{ID} = $parentid_1 * $parentid_2;
        $child{ParentID} = $parentid_1 . "_" . $parentid_2;
        $child{age} = 0;
        $child{status} = 0;
        $child{matepair} = 0;
        $child{matenum} = 0;

        #子供の遺伝子型生成
        #XY
        my $XY_1 = $$in{$id}{XYtype};
        my $XY_2 = $$in{$pair_id}{XYtype};
        my $gamete_xy_1 = substr($XY_1, int(rand(1)+0.5), 1);
        my $gamete_xy_2 = substr($XY_2, int(rand(1)+0.5), 1);
        my @gamete_xy = ($gamete_xy_1, $gamete_xy_2);
        @gamete_xy = sort{$a cmp $b} @gamete_xy;
        $child{XYtype} = join("",@gamete_xy);
        my %code_xy = ("XX" => 1, "XY" => 0, "YY" => -1);
        $child{XYcode} = $code_xy{$child{XYtype}};

        #ABO
        my $ABO_1 = $$in{$id}{ABO};
        my $ABO_2 = $$in{$pair_id}{ABO};
        my $gamete_abo_1 = substr($ABO_1, int(rand(1)+0.5), 1);
        my $gamete_abo_2 = substr($ABO_2, int(rand(1)+0.5), 1);
        my @gamete_abo = ($gamete_abo_1, $gamete_abo_2);
        @gamete_abo = sort{$a cmp $b} @gamete_abo;
        $child{ABO} = join("", @gamete_abo);
        my $omegatype = &Omegatype(\%child);
        $child{omegatype} = $omegatype;
        my %code_abo = ("Alpha" => 2, "Beta" => 1, "Omega" => 0);
        $child{omegacode} = $code_abo{$omegatype};

        #個体群に代入
        $population++;
        $Omegaverse->{$population} = \%child;
    }
}

#死亡判定・Omegaverseの通し番号整理・Betaのつがい処理
sub Death {
    my $out;
    my $in = $_[0];
    my @pop = sort keys(%$in);
    my $new_id = 0;
    my %hash;
    foreach my $id (sort keys($in)){
        my %Individual;

        #死亡処理：男女別に死力表読み込み生成乱数から処理
        my $age= $$in{$id}{age} + 1;
        my $die = rand(1);
        if($$in{$id}{XYcode} eq 1){
            next if $die < $death_XX{$age};
        }elsif($$in{$id}{XYcode} eq 0){
            next if $die < $death_XY{$age};
        #YY個体は致死
        }else{
            next;
        }
        #各種項目をコピー
        foreach my $key (keys(%$in{$id})){
            $Individual{$key} = $$in{$id}{$key};
        }
        #交配可能判定(status:2はここでリセット)
        my $status = &Matestatus(\%Individual);
        $Individual{status} = $status;

        #新規個体群表に代入
        $new_id++;
        $hash{$id} = $new_id;
        $out->{$new_id} = \%Individual;
    }

    #matenumの旧idから新idの付け替え(たぶんもっと上手い方法がある)
    foreach my $id (sort keys($out)){
        my $old_mate = $$out{$id}{matenum};
        #つがい死亡時の処理(betaなら0に、それ以外は負数に(Matestatus用の処理))
        next if $old_mate eq 0;
        unless(defined($hash{$old_mate})){
            if($$out{$id}{omegatype} eq "Beta"){
                $$out{$id}{matenum} = 0;
                $$out{$id}{matepair} = 0;
                my $status = &Matestatus(%$out{$id});
                $$out{$id}{status} = $status;
            }
        }else{
            $$out{$id}{matenum} = $hash{$old_mate};
        }
    }
    return($out);
}
#####################################
